import os, time, json, shutil, subprocess
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS

def _has_wcs(header) -> bool:
    try:
        w = WCS(header); return w.has_celestial
    except Exception:
        return False

def _fits_with_wcs(path: Path) -> bool:
    with fits.open(path) as hdul:
        for h in hdul:
            if h.data is not None and getattr(h.data, 'ndim', 0) == 2:
                if _has_wcs(h.header): return True
    return False

def ensure_wcs(input_path: Path, out_dir: Path) -> Path:
    import requests
    def log(s): print(f"[astro-sky-solver] {s}")
    if input_path.suffix.lower() in ('.fits','.fit','.fts') and _fits_with_wcs(input_path):
        return input_path
    api_key = os.environ.get('ASTROMETRY_API_KEY','').strip()
    if api_key:
        log('Using astrometry.net web API...')
        r = requests.post('https://nova.astrometry.net/api/login', data={'request-json': json.dumps({'apikey': api_key})}, timeout=30)
        payload = r.json()
        if payload.get('status') == 'error' or 'session' not in payload:
            raise RuntimeError(f'API login error: {payload}')
        sess = payload['session']
        log('Uploading image...')
        with open(input_path,'rb') as fh:
            files={'file': fh}
            r = requests.post('https://nova.astrometry.net/api/upload', data={'request-json': json.dumps({'session': sess})}, files=files, timeout=180)
        up = r.json()
        if up.get('status') == 'error' or not up.get('subid'): raise RuntimeError(f'Upload error: {up}')
        subid = up['subid']
        job_id=None
        for _ in range(120):
            time.sleep(5)
            rr = requests.get(f'https://nova.astrometry.net/api/submissions/{subid}', timeout=30)
            jobs = [j for j in rr.json().get('jobs',[]) if j is not None]
            if jobs: job_id=jobs[0]; break
            log('Waiting for job assignment...')
        if not job_id: raise RuntimeError('Timed out waiting for job id.')
        while True:
            time.sleep(5)
            jr = requests.get(f'https://nova.astrometry.net/api/jobs/{job_id}', timeout=30)
            status = jr.json().get('status','')
            log(f'Job {job_id} status: {status}')
            if status in ('success','failure'): break
        if status != 'success': raise RuntimeError('Plate solve failed on server.')
        solved_url = f'https://nova.astrometry.net/new_fits_file/{job_id}'
        out = out_dir / (input_path.stem + '_wcs.fits')
        resp = requests.get(solved_url, timeout=180, allow_redirects=True)
        if not resp.ok or not resp.content.startswith(b'SIMPLE'): raise RuntimeError('Failed to download solved FITS.')
        with open(out,'wb') as f: f.write(resp.content)
        log(f'Downloaded solved FITS with WCS: {solved_url}')
        return out
    if shutil.which('solve-field'):
        log('Using local astrometry.net solve-field...')
        out_base = out_dir / (input_path.stem + '_wcs')
        cmd = ['solve-field', str(input_path),'--overwrite','--no-plots','--downsample','2','--dir',str(out_dir),'--new-fits',str(out_base.with_suffix('.fits'))]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode!=0: print(res.stdout, res.stderr); raise RuntimeError('solve-field failed.')
        out = out_base.with_suffix('.fits')
        if not out.exists(): raise RuntimeError('Expected WCS FITS missing.')
        return out
    raise RuntimeError('Need WCS: provide solved FITS, set ASTROMETRY_API_KEY, or install solve-field.')
