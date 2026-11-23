#!/usr/bin/env python3
"""Prepare red/DM noise (without an early nuisance GP), sample/reuse and verify."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

# Set BLAS limits before subprocesses import NumPy (including multiprocessing).
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
SCRIPTS=Path(__file__).resolve().parent/'scripts'
sys.path.insert(0,str(SCRIPTS))
from common import read_config, json_write, sha256, snapshot_tim


def nonnegative_int(value):
    number=int(value)
    if number<0:
        raise argparse.ArgumentTypeError('Require a nonnegative integer (0 disables the GP)')
    return number


def validate_output_directory(out, protected_paths):
    """Reject output paths that could overwrite inputs or follow a symlink."""
    if out.is_symlink():
        raise ValueError('outdir must not be a symbolic link')
    out=out.resolve()
    for protected in [Path.cwd(),Path.home(),Path(__file__).resolve().parent,
                      *protected_paths]:
        protected=Path(protected).resolve()
        if out==protected or out in protected.parents:
            raise ValueError(f'outdir would overwrite an input or working directory: {protected}')
    if out.exists() and not out.is_dir():
        raise ValueError('outdir exists but is not a directory')


def remove_output(path):
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    elif path.exists() or path.is_symlink():
        path.unlink()


def reset_output_directory(out, protected_paths):
    """Start a fresh run without deleting inputs or following directory symlinks."""
    validate_output_directory(out,protected_paths)
    out=out.resolve()
    if out.exists():
        print('Overwriting output directory:',out,flush=True)
        for child in out.iterdir():
            remove_output(child)
    else:
        out.mkdir(parents=True)


def export_results(results, out):
    """Expose the final result files and figures directly in the requested outdir."""
    for source in results.iterdir():
        destination=out/source.name
        if source.is_dir():
            shutil.copytree(source,destination,dirs_exist_ok=True)
        else:
            shutil.copy2(source,destination)
    json_write(out/'exported_products.json',[source.name for source in results.iterdir()])


def clear_exported_results(out):
    """Remove obsolete top-level results while retaining the prepared data/chain."""
    products=out/'exported_products.json'
    names=json.loads(products.read_text()) if products.exists() else []
    final=out/'final_result.json'
    if not names and final.exists():
        previous=json.loads(final.read_text())
        source=Path(previous.get('pass_results_directory',previous['results_directory']))
        if source.resolve()!=out and source.is_dir():
            names=[item.name for item in source.iterdir()]
    for name in names:
        if Path(name).name!=name or name in ('.','..','logs') or name.startswith('pass_'):
            raise ValueError('Invalid exported product name')
    for name in names:
        remove_output(out/name)
    remove_output(final)
    remove_output(products)


def latest_chain_pass(out):
    candidates=[]
    if out.is_dir():
        for folder in out.iterdir():
            match=re.fullmatch(r'pass_(\d+)',folder.name)
            if match and (folder/'noise'/'chain.h5').exists():
                candidates.append(int(match[1]))
    return max(candidates,default=None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('parfile',type=Path)
    parser.add_argument('timfile',type=Path)
    parser.add_argument('--outdir',required=True,type=Path,
                        help='Output directory; automatically resume an existing MCMC chain, or overwrite when starting fresh/reusing samples')
    parser.add_argument('--restart',action='store_true',
                        help='Clear existing output (including all MCMC chains) and start again')
    parser.add_argument('--config',type=Path)
    parser.add_argument('--dm-components','--dm-modes',type=nonnegative_int,
                        help='DM Fourier mode count; overrides config dm_modes (default: config value or 20; 0 disables DM GP)')
    parser.add_argument('--red-components','--red-modes',type=nonnegative_int,
                        help='Red Fourier mode count; overrides config red_modes (default: config value or 20; 0 disables red GP)')
    custom=parser.add_mutually_exclusive_group()
    custom.add_argument('--use-ddgrab',action='store_true')
    custom.add_argument('--use-ddab',action='store_true')
    custom.add_argument('--use-ddhab',action='store_true')
    shapiro=parser.add_mutually_exclusive_group()
    shapiro.add_argument('--sample-dd-shapiro',action='store_true')
    shapiro.add_argument('--sample-ddh-shapiro',action='store_true')
    parser.add_argument('--reuse-noise',type=Path,help='Directory with posterior_samples.npy and pars.txt')
    parser.add_argument('--reuse-mode',choices=['auto','identical','reweight','conditional'],default='auto')
    parser.add_argument('--old-data',type=Path,help='Source likelihood_data.npz for importance reweighting')
    parser.add_argument('--initial-samples',type=Path,help='Initialize a fresh MCMC from paired old draws without changing priors')
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--moves',choices=['mixed','mixed-prior','mixed-white','stretch'],default='mixed')
    parser.add_argument('--walkers',type=int,default=48)
    parser.add_argument('--tau-multiple',type=float,default=50.)
    parser.add_argument('--tau-stability',type=float,default=.02)
    parser.add_argument('--min-steps',type=int,default=5000,
                        help='Minimum total steps per chain, including steps already saved')
    parser.add_argument('--max-steps',type=int,default=30000,
                        help='Maximum total steps per chain; increase this cap to extend an existing chain')
    parser.add_argument('--draws',type=int,default=4000)
    parser.add_argument('--max-passes',type=int,default=3)
    parser.add_argument('--allow-downloads',action='store_true')
    args=parser.parse_args()
    if args.max_passes<1 or args.draws<20:
        raise ValueError('Require at least one pass and 20 recovered draws')
    config=read_config(args.config)
    if args.dm_components is not None: config['dm_modes']=args.dm_components
    if args.red_components is not None: config['red_modes']=args.red_components
    identity=dict(par_sha256=sha256(args.parfile),tim_sha256=sha256(args.timfile),
                  config=config,use_ddgrab=args.use_ddgrab,sample_dd_shapiro=args.sample_dd_shapiro,
                  reuse_noise=str(args.reuse_noise.resolve()) if args.reuse_noise else None,
                  reuse_mode=args.reuse_mode,
                  reuse_samples_sha256=sha256(args.reuse_noise/'posterior_samples.npy') if args.reuse_noise else None,
                  old_data=str(args.old_data.resolve()) if args.old_data else None)
    if args.use_ddab: identity['use_ddab']=True
    if args.use_ddhab: identity['use_ddhab']=True
    if args.sample_ddh_shapiro: identity['sample_ddh_shapiro']=True
    # Resolve all input dependencies before removing an earlier run, including
    # TIM INCLUDE files and the old likelihood referenced by supplied samples.
    protected=[args.parfile,args.timfile]
    for path in (args.config,args.reuse_noise,args.old_data,args.initial_samples):
        if path is not None: protected.append(path)
    with tempfile.TemporaryDirectory(prefix='binary-noise-input-check-') as temporary:
        manifest=snapshot_tim(args.timfile,Path(temporary)/'input.tim')
        protected.extend(entry['original'] for entry in manifest)
    for source in (args.reuse_noise,args.initial_samples):
        if source is not None and (source/'source_data.json').exists():
            source_data=json.loads((source/'source_data.json').read_text())['datafile']
            protected.append(Path(source_data).resolve().parent)
    if args.old_data: protected.append(args.old_data.resolve().parent)
    validate_output_directory(args.outdir,protected)
    out=args.outdir.resolve()
    signature=out/'run_identity.json'
    chain_pass=latest_chain_pass(out) if not args.reuse_noise and not args.restart else None
    resuming=chain_pass is not None
    if resuming:
        if not signature.exists() or json.loads(signature.read_text())!=identity:
            raise ValueError('Existing MCMC inputs/configuration differ or lack provenance; use --restart to clear the old chain')
        data=out/f'pass_{chain_pass:02d}'/'prepared'
        noise=out/f'pass_{chain_pass:02d}'/'noise'
        if any(not path.resolve().is_relative_to(out) for path in (data,noise,noise/'chain.h5')):
            raise ValueError('Existing chain/prepared paths must remain inside outdir')
        if not (data/'likelihood_data.npz').exists() or not (data/'metadata.json').exists():
            raise ValueError('Existing chain lacks prepared data/metadata; use --restart')
        metadata=json.loads((data/'metadata.json').read_text())
        if args.use_ddgrab and not metadata.get('ddgrab_fit_flags_preserved',False):
            raise ValueError('Existing DDGRAB chain used the legacy forced fit flags; use --restart to prepare data with the input par flags')
        if metadata.get('allow_downloads',False)!=args.allow_downloads:
            raise ValueError('--allow-downloads differs from the existing run; keep it unchanged or use --restart')
        for archived in metadata['tim_manifest']:
            if sha256(archived['original'])!=archived['original_sha256']:
                raise ValueError('A TIM INCLUDE changed; use --restart to clear the old chain')
        for key,value in metadata.get('runtime_overrides',{}).items():
            if os.environ.get(key)!=value:
                raise ValueError(f'{key} differs from the existing chain; use the original environment or --restart')
        previous_arguments=out/'command_arguments.json'
        if previous_arguments.exists() and json.loads(previous_arguments.read_text())['walkers']!=args.walkers:
            raise ValueError('Walker count differs from the existing chain; keep --walkers unchanged or use --restart')
        if args.max_passes<chain_pass:
            raise ValueError(f'--max-passes must be at least {chain_pass} to resume this chain')
        print('Resuming MCMC chain:',out/f'pass_{chain_pass:02d}'/'noise'/'chain.h5',flush=True)
        clear_exported_results(out)
        # Later passes without a chain can contain incomplete/stale preparation.
        for folder in out.iterdir():
            match=re.fullmatch(r'pass_(\d+)',folder.name)
            if match and int(match[1])>chain_pass: remove_output(folder)
    else:
        reset_output_directory(args.outdir,protected)
    logs=out/'logs'; logs.mkdir(exist_ok=True)
    json_write(signature,identity)
    json_write(out/'resolved_config.json',config)
    json_write(out/'command_arguments.json',{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()})
    def run(script,arguments,label):
        command=[sys.executable,'-B',str(SCRIPTS/script)]+list(map(str,arguments))
        print('Running:', ' '.join(command),flush=True)
        with (logs/(label+'.log')).open('a' if resuming else 'w') as stream:
            result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f'{label} failed (exit {result.returncode}); see {logs/(label+".log")}')
    options=[]
    if args.use_ddgrab: options+=['--use-ddgrab']
    if args.use_ddab: options+=['--use-ddab']
    if args.use_ddhab: options+=['--use-ddhab']
    if args.sample_dd_shapiro: options+=['--sample-dd-shapiro']
    if args.sample_ddh_shapiro: options+=['--sample-ddh-shapiro']
    if args.allow_downloads: options+=['--allow-downloads']
    par=args.parfile.resolve()
    first_pass=chain_pass if resuming else 1
    passes_file=out/'passes.json'
    passes=[p for p in json.loads(passes_file.read_text()) if p['pass_number']<first_pass] if resuming and passes_file.exists() else []
    for number in range(first_pass,args.max_passes+1):
        folder=out/f'pass_{number:02d}'
        data=folder/'prepared'; noise=folder/'noise'; results=folder/'results'
        if not (resuming and number==first_pass):
            run('prepare.py',[par,args.timfile.resolve(),'--outdir',data,'--config',out/'resolved_config.json']+options,f'prepare_{number:02d}')
        if args.reuse_noise:
            arguments=[data/'likelihood_data.npz',args.reuse_noise.resolve(),'--outdir',noise,'--mode',args.reuse_mode]
            if args.old_data: arguments+=['--old-data',args.old_data.resolve()]
            run('reuse.py',arguments,f'reuse_{number:02d}')
        else:
            run('sample.py',[data/'likelihood_data.npz','--outdir',noise,'--workers',args.workers,
                 '--walkers',args.walkers,'--moves',args.moves,'--tau-multiple',args.tau_multiple,'--tau-stability',args.tau_stability,'--min-steps',args.min_steps,'--max-steps',args.max_steps,
                 '--seed',161916+number]+(['--initial-samples',args.initial_samples.resolve()] if args.initial_samples else []),f'sample_{number:02d}')
        remove_output(results)
        run('recover.py',[data,noise,'--outdir',results,'--draws',args.draws],f'recover_{number:02d}')
        summary=json.loads((results/'analysis_summary.json').read_text())
        passes.append(dict(pass_number=number,results=str(results),accepted=summary['numerically_accepted']))
        json_write(out/'passes.json',passes)
        if summary['numerically_accepted']:
            export_results(results,out)
            json_write(out/'final_result.json',dict(**summary,results_directory=str(out),
                                                  pass_results_directory=str(results)))
            print('Completed:',out,flush=True)
            return
        par=results/'result.par'
    raise RuntimeError('Local timing approximation did not pass verification; all passes retained. Inspect parameterization/model before accepting.')


if __name__=='__main__':
    main()
