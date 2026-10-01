"""Publication plots for the joint fit; native per-setting plots plus joint QA.

Called by run_joint_xsec_fit.py. ROOT is loaded before building the native
extractor in render-only mode. Additional figures use Matplotlib's Agg backend.
"""
import json
import math
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np

import run_joint_xsec_fit as fit


def run_native(config, vertex, folder, env, partons=False, warmups=10000, calls=100000,
               snapshot=None, positive=True, variance='data'):
    """Build the native extractor privately; the reference shell stays read-only."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    flags = shlex.split(subprocess.check_output(['root-config', '--cflags', '--libs'], env=env, text=True))
    with tempfile.TemporaryDirectory(prefix='joint_plot_build_') as temporary:
        build = Path(temporary)
        (build/'xsec_config.h').write_text(fit.render(config))
        binary = build/'joint_plot_renderer'
        command = [env.get('CXX', 'g++'), '-std=c++17', '-O2', '-I'+str(build), '-I'+str(fit.HERE),
                   str(fit.HERE/'excl_xsec_pi0_analysis_no_simc_model.C'), *flags, '-lMinuit2']
        native_env = dict(env)
        if partons:
            software = Path(env.get('NPS_SOFTWARE_ROOT', '/group/nps/singhav/software'))
            model = Path(env.get('NPS_PARTONS_ROOT', str(software/'partons')))
            libs = [model/'lib64/libsfml-system.so', model/'lib/libcln.so',
                    model/'lib/libElementaryUtils.so', model/'lib/libNumA++.so',
                    model/'lib/libPARTONS.so', software/'python/lib/libgsl.so',
                    software/'python/lib/libgslcblas.so', software/'apfel/lib/libapfelxx.so',
                    software/'lhapdf/lib/libLHAPDF.so', Path('/usr/lib64/libxml2.so'),
                    software/'lhapdf/lib/libstdc++.so']
            schema = software/'src/partons/partons-example/data/xmlSchema.xsd'
            for path in [*libs, schema]:
                fit.need(path.is_file(), f'missing native PARTONS dependency: {path}')
            (build/'partons.properties').write_text(
                f'log.file.path = {build}/logger.properties\nxml.schema.file.path = {schema}\n'
                'computation.nb.processor = 1\ngpd.service.batch.size = 1000\n'
                'collinear_distribution.service.batch.size = 1000\nccf.service.batch.size = 1000\n'
                'observable.service.batch.size = 1000\n')
            (build/'logger.properties').write_text(
                f'enable = true\ndefault.level = WARN\nprint.mode = COUT\nlog.folder.path = {build}\n')
            libdirs = [model/'lib64', model/'lib', software/'python/lib', software/'apfel/lib', software/'lhapdf/lib']
            rpath = ':'.join(str(p) for p in libdirs)
            command += ['-DNPS_ENABLE_PARTONS', '-I'+str(model/'include'),
                        '-I'+str(software/'python/include'), '-I'+str(software/'apfel/include'),
                        '-I'+str(software/'lhapdf/include'), '-I/usr/include/libxml2',
                        '-Wl,-rpath,'+rpath, *(str(p) for p in libs)]
            native_env['LD_LIBRARY_PATH'] = str(software/'lhapdf/lib')+':'+rpath+':'+env.get('LD_LIBRARY_PATH', '')
        command += ['-o', str(binary)]
        run = [str(binary), '--data-file', config['data_file'], '--sim-file', config['simc_file'],
               '--vertex_simc_file', str(vertex), '--kin', config['configured_kinematic'],
               '--out-dir', str(folder), '--out-root', str(folder/'joint_plot_diagnostics.root'),
               '--out-slice-csv', str(folder/'joint_plot_slices.csv'),
               '--all-plots-pdf', str(folder/'all_setting_plots.pdf'),
               '--fit-objective', 'gaussian', '--fit-variance', variance,
               '--positive-xsec' if positive else '--no-positive-xsec']
        if snapshot is not None:
            run += ['--joint-plot-input', str(snapshot)]
        if config['mmiss_select'] != 'window':
            run += ['--mmiss_select', config['mmiss_select']]
        if partons:
            run += ['--partons', '--partons-warmups', str(warmups), '--partons-calls', str(calls)]
        with (folder/'render.log').open('w') as log:
            try:
                subprocess.run(command, env=native_env, cwd=fit.REPO, stdout=log, stderr=subprocess.STDOUT, check=True)
                subprocess.run(run, env=native_env, cwd=fit.REPO, stdout=log, stderr=subprocess.STDOUT, check=True)
            except subprocess.CalledProcessError as error:
                raise ValueError(f"native plot rendering failed; see {folder/'render.log'}") from error


def matrix(path, size):
    result = np.full((size, size), np.nan)
    for row in fit.rows(path):
        result[int(row['parameter_i']), int(row['parameter_j'])] = float(row['stat_plus_mc_covariance'])
    return result


def setting_snapshot(output, index, setting, destination):
    """Local marginal covariance is sufficient for that setting's fit curves."""
    parameters = fit.rows(output/'joint_parameters.csv')
    keys = {(int(r['truth_block']), int(r['setting_index']), r['component']): i
            for i, r in enumerate(parameters)}
    covariance = matrix(output/'joint_covariance.csv', len(parameters))
    blocks = sorted(b for b, s, c in keys if s == index and c == 'U')
    columns = [keys[b, index if c == 'U' else -1, c] for b in blocks for c in fit.COMPONENTS]
    joint_rows = [r for r in fit.rows(output/'joint_rows.csv') if int(r['setting_index']) == index]
    joint_rows.sort(key=lambda r: int(r['reco_row']))
    cells, _, _ = fit.load_response_cells(Path(setting['output_dir']))
    cells = {(int(r['reco_row']), int(r['truth_block'])): r for r in cells}
    summary = dict(line.split('=', 1) for line in (output/'joint_summary.txt').read_text().splitlines())
    positivity = {int(r['truth_block']): r for r in fit.rows(output/'joint_positivity.csv')
                  if int(r['setting_index']) == index}
    with destination.open('w') as stream:
        def line(values):
            stream.write(' '.join(str(v) for v in values)+'\n')
        line(['joint_plot_v1', len(blocks), len(joint_rows), summary['positivity_boundary_active']])
        line([summary[k] for k in ('chi2', 'ndf', 'rank', 'condition', 'mc_iterations')])
        for b in blocks:
            tolerance = float(positivity[b]['feasibility_tolerance'])
            line([b, tolerance if math.isfinite(tolerance) else 0.] +
                 [parameters[keys[b, index if c == 'U' else -1, c]]['value'] for c in fit.COMPONENTS])
        for i in columns:
            line(covariance[i, columns])
        for r, row in enumerate(joint_rows):
            fit.need(int(row['reco_row']) == r, 'joint plotting rows must be contiguous')
            line([int(int(row['fit_index']) >= 0), row['data'], row['data_variance'],
                  row['variance_used'], row['prediction']] +
                 [cells[r, b][name] for b in blocks for name in fit.BASIS])


def native_plots(output, manifest, partons, warmups, calls):
    env = fit.root_environment()
    # The saved fit, not current extraction environment overrides, owns inputs.
    env = {k: v for k, v in env.items() if not k.startswith('NPS_XSEC_')}
    pages = []
    for index, setting in enumerate(manifest['settings']):
        folder = output/'plots'/'settings'/setting['kinematic']
        folder.mkdir(parents=True, exist_ok=True)
        snapshot = folder/'joint_plot_input.txt'
        setting_snapshot(output, index, setting, snapshot)
        config_path, config, _ = fit.read_config(Path(setting['config']))
        source = Path(setting['output_dir'])
        meta_path = source/'joint_input_metadata.txt'
        metadata = fit.read_metadata(meta_path if meta_path.is_file() else source/Path(config['out_root']).name)
        vertex = metadata.get('vertex_source')
        fit.need(vertex and Path(vertex).is_file(),
                 f"{setting['kinematic']}: original vertex_source is required to rebuild all pipeline plots")
        # Saved extraction cuts can differ from a preset's defaults. Freeze the
        # renderer to the actual prepared/fitted selections and normalization.
        config.update(mmiss_select=manifest['selection']['exclusive_selection'],
                      mmiss_lower_gev=float(manifest['selection']['mmiss_lower_gev']),
                      mmiss_upper_gev=float(manifest['selection']['mmiss_upper_gev']),
                      tgt_contam=float(manifest['selection']['target_contam_factor']),
                      tgt_contam_err=float(manifest['selection']['target_contam_factor_err']),
                      data_file=setting['input_data_file'], simc_file=setting['input_simc_file'])
        effective = folder/'plot_config.json'
        effective.write_text(json.dumps(config, indent=2)+'\n')
        print(f"[PLOTS] {setting['kinematic']}: native pipeline plots from joint results", flush=True)
        run_native(config, vertex, folder, env, partons, warmups, calls, snapshot,
                   manifest['options']['positive_xsec'], manifest['options']['fit_variance'])
        native_pages = sorted(p for p in folder.rglob('*.pdf') if p.name != 'all_setting_plots.pdf')
        fit.need(native_pages and (folder/'joint_plot_slices.csv').is_file() and
                 (folder/'all_setting_plots.pdf').is_file(),
                 f"{setting['kinematic']}: native plot rendering produced incomplete outputs")
        pages.extend(native_pages)
    return pages


def joint_figures(output, manifest):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    plt.rcParams.update({'font.size': 11, 'axes.grid': True, 'grid.alpha': .2,
                         'savefig.dpi': 160, 'figure.constrained_layout.use': True})
    folder = output/'plots'/'joint'
    folder.mkdir(parents=True, exist_ok=True)
    pages = []
    scale = 1.e9  # Native SIMC microbarn/MeV^2 -> nb/GeV^2, as in native plots.
    unit = r'nb/GeV$^2$'
    parameters = fit.rows(output/'joint_parameters.csv')
    separated = fit.rows(output/'joint_separated_parameters.csv')
    cov = matrix(output/'joint_covariance.csv', len(parameters))
    sep_cov = matrix(output/'joint_separated_covariance.csv', len(separated))
    keys = {(int(r['truth_block']), int(r['setting_index']), r['component']): i
            for i, r in enumerate(parameters)}
    sepkeys = {(int(r['truth_block']), r['component']): i for i, r in enumerate(separated)}
    nominal = manifest['lt_separation']['settings']
    names = [s['kinematic'] for s in nominal]
    blocks = sorted({int(r['truth_block']) for r in separated if r['region'] == 'published'})
    groups = sorted({(int(r['iq']), int(r['ix'])) for r in separated if r['region'] == 'published'})
    bins = manifest['bins']
    statuses = {r['truth_block']: r['status'] for r in manifest['lt_separation']['blocks']}

    def save(fig, name):
        path = folder/(name+'.pdf')
        fig.savefig(path)
        fig.savefig(folder/(name+'.png'))
        pages.append(path)
        plt.close(fig)

    def values(records, indices):
        return np.array([float(records[i]['value']) for i in indices])*scale

    def errorbar(ax, x, y, variance, **kwargs):
        # Unavailable errors remain visibly central-only, never zero errors.
        x, y, variance = map(np.asarray, (x, y, variance))
        valid = np.isfinite(x) & np.isfinite(y)
        errors = valid & np.isfinite(variance) & (variance >= 0)
        color = kwargs.pop('color', None)
        label = kwargs.pop('label', None)
        if errors.any():
            artist = ax.errorbar(x[errors], y[errors], yerr=np.sqrt(variance[errors]),
                                fmt='o', capsize=3, label=label, color=color, **kwargs)
            color = artist[0].get_color()
        missing = valid & ~errors
        if missing.any():
            ax.plot(x[missing], y[missing], 'x', color=color,
                    label=(label+'; ' if label else '')+'central only', **kwargs)

    def block_label(b):
        r = separated[sepkeys[b, 'T']]
        it, iq, ix = (int(r[k]) for k in ('it', 'iq', 'ix'))
        return (f"Truth bin {b}: " +
                rf"$-t'\in[{0. - bins['tprime_bin_edges'][it+1]:g},{0. - bins['tprime_bin_edges'][it]:g}]$ GeV$^2$" +
                '\n'+rf"$Q^2\in[{bins['q2_bin_edges'][iq]:g},{bins['q2_bin_edges'][iq+1]:g}]$, " +
                rf"$x_B\in[{bins['xb_bin_edges_by_q2'][iq][ix]:g},{bins['xb_bin_edges_by_q2'][iq][ix+1]:g}]$")

    for b in blocks:
        fig, ax = plt.subplots(figsize=(8, 6))
        uidx = [keys[b, s, 'U'] for s in range(len(names)) if (b, s, 'U') in keys]
        settings = [int(parameters[i]['setting_index']) for i in uidx]
        eps = np.array([nominal[s]['epsilon'] for s in settings])
        for i, s, e in zip(uidx, settings, eps):
            errorbar(ax, [e], values(parameters, [i]), [cov[i, i]*scale**2], label=names[s])
        ti, li = sepkeys[b, 'T'], sepkeys[b, 'L']
        t, longitudinal = values(separated, [ti, li])
        if np.isfinite([t, longitudinal]).all():
            grid = np.linspace(0, 1, 200)
            ax.plot(grid, t+grid*longitudinal, color='black', label=r'$\sigma_T+\epsilon\sigma_L$')
            c = sep_cov[np.ix_([ti, li], [ti, li])]*scale**2
            if np.isfinite(c).all():
                design = np.column_stack((np.ones(len(grid)), grid))
                err = np.sqrt(np.maximum(0, np.einsum('ij,jk,ik->i', design, c, design)))
                ax.fill_between(grid, t+grid*longitudinal-err, t+grid*longitudinal+err,
                                color='gray', alpha=.2, label='Marginal 1 SD band')
        ax.set(xlabel=r'Nominal $\epsilon$', ylabel=rf'$\sigma_U$ [{unit}]',
               title=block_label(b)+'\n'+statuses[b], xlim=(0, 1))
        ax.legend(fontsize=9)
        save(fig, f'rosenbluth_truth_{b}')

        fig, axes = plt.subplots(1, 2, figsize=(11, 5))
        for ax, components in zip(axes, [('T', 'L'), ('LT', 'TT')]):
            ii = [sepkeys[b, c] for c in components]
            center = values(separated, ii)
            subcov = sep_cov[np.ix_(ii, ii)]*scale**2
            if np.isfinite(center).all():
                ax.plot(*center, 'o', color='black')
                if np.isfinite(subcov).all():
                    eigen, vectors = np.linalg.eigh(subcov)
                    if (eigen >= 0).all():
                        angle = np.degrees(np.arctan2(vectors[1, 1], vectors[0, 1]))
                        radius = np.sqrt(2.30*eigen)
                        ax.add_patch(Ellipse(center, 2*radius[1], 2*radius[0], angle=angle,
                                             facecolor='C0', edgecolor='C0', alpha=.25))
                        ax.update_datalim([center-1.8*np.sqrt(np.diag(subcov)),
                                           center+1.8*np.sqrt(np.diag(subcov))])
                        ax.autoscale_view()
                else:
                    ax.text(.5, .9, 'Covariance unavailable; central only', transform=ax.transAxes, ha='center', fontsize=9)
            else:
                ax.text(.5, .5, statuses[b], transform=ax.transAxes, ha='center')
            ax.set(xlabel=rf'$\sigma_{{{components[0]}}}$ [{unit}]', ylabel=rf'$\sigma_{{{components[1]}}}$ [{unit}]')
        fig.suptitle(block_label(b)+'\nGaussian 68% joint ellipses (two parameters; conditional covariance)')
        save(fig, f'coefficient_ellipses_truth_{b}')

    for iq, ix in groups:
        chosen = [b for b in blocks if int(separated[sepkeys[b, 'T']]['iq']) == iq
                  and int(separated[sepkeys[b, 'T']]['ix']) == ix]
        chosen.sort(key=lambda b: int(separated[sepkeys[b, 'T']]['it']))
        x = [-.5*(bins['tprime_bin_edges'][int(separated[sepkeys[b, 'T']]['it'])] +
                   bins['tprime_bin_edges'][int(separated[sepkeys[b, 'T']]['it'])+1]) for b in chosen]
        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        for ax, component in zip(axes.flat, ('T', 'L', 'LT', 'TT')):
            ii = [sepkeys[b, component] for b in chosen]
            errorbar(ax, x, values(separated, ii), np.diag(sep_cov)[ii]*scale**2)
            ax.axhline(0, color='gray', lw=.8)
            ax.set(xlabel=r"$-t'$ [GeV$^2$] (bin centers)", ylabel=rf'$\sigma_{{{component}}}$ [{unit}]')
            if not np.isfinite(values(separated, ii)).all():
                ax.text(.03, .96, 'Unavailable bins omitted; see separation status', va='top', transform=ax.transAxes, fontsize=8)
        fig.suptitle(f'Joint coefficients and post-fit L/T separation: Q2 bin {iq}, xB bin {ix}\nMarginal stat + MC errors; fixed nominal epsilon')
        save(fig, f'separated_terms_q{iq}_x{ix}')
        fig, ax = plt.subplots(figsize=(8, 6))
        for s, name in enumerate(names):
            ii = [keys[b, s, 'U'] for b in chosen]
            errorbar(ax, x, values(parameters, ii), np.diag(cov)[ii]*scale**2, label=name)
        ax.set(xlabel=r"$-t'$ [GeV$^2$] (bin centers)", ylabel=rf'$\sigma_U$ [{unit}]', title=f'Independent U values: Q2 bin {iq}, xB bin {ix}')
        ax.legend(fontsize=9)
        save(fig, f'per_setting_U_q{iq}_x{ix}')

    for title, records, covariance in [('joint', parameters, cov), ('separated', separated, sep_cov)]:
        sd = np.sqrt(np.maximum(0, np.diag(covariance)))
        with np.errstate(divide='ignore', invalid='ignore'):
            correlation = covariance / np.outer(sd, sd)
        fig, ax = plt.subplots(figsize=(10, 9))
        palette = plt.get_cmap('RdBu_r').copy()
        palette.set_bad('#dedede')
        im = ax.imshow(np.ma.masked_invalid(correlation), vmin=-1, vmax=1, cmap=palette)
        labels = [f"b{r['truth_block']}:{r['component']}" +
                  (f":s{r['setting_index']}" if r.get('component') == 'U' else '') for r in records]
        step = max(1, math.ceil(len(labels)/40))
        ticks = list(range(0, len(labels), step))
        ax.set_xticks(ticks, [labels[i] for i in ticks], rotation=90, fontsize=7)
        ax.set_yticks(ticks, [labels[i] for i in ticks], fontsize=7)
        ax.set_title(f'{title.capitalize()} parameter correlations\nGray = covariance unavailable; guards included')
        fig.colorbar(im, ax=ax, label='Correlation')
        save(fig, title+'_correlation')

    records = fit.rows(output/'joint_rows.csv')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for s, name in enumerate(names):
        selected = [r for r in records if int(r['setting_index']) == s and int(r['fit_index']) >= 0]
        r = np.array([int(v['reco_row']) for v in selected])
        pulls = np.array([float(v['pull']) for v in selected])
        axes[0, 0].plot(r, pulls, '.', label=name)
        axes[0, 1].hist(pulls, bins=np.linspace(-5, 5, 31), histtype='step', label=name)
        fraction = [float(v['mc_variance_final'])/float(v['variance_used']) for v in selected]
        axes[1, 0].plot(r, fraction, '.', label=name)
        axes[1, 1].bar(s, np.sum(pulls**2), label=name)
    axes[0, 0].axhline(0, color='gray')
    axes[0, 0].set(xlabel='Reconstructed row', ylabel=r'$(data-prediction)/\sqrt{V_{fit}}$')
    axes[0, 1].set(xlabel='Fit-standardized residual (not leverage corrected)', ylabel='Rows in [-5,5]')
    axes[1, 0].set(xlabel='Reconstructed row', ylabel='MC variance / fit variance',
                   title='MC variance diagnostic' if manifest['options']['fit_variance'] == 'finite-mc' else 'MC variance shown but excluded from data-only fit')
    axes[1, 1].set(xticks=range(len(names)), xticklabels=[f's{s}' for s in range(len(names))], ylabel=r'$\chi^2$ contribution')
    axes[0, 0].legend(fontsize=8)
    fig.suptitle('Joint fit residuals and per-setting contributions\nResiduals are correlated through the fit; no per-setting degrees of freedom assigned')
    save(fig, 'joint_fit_quality')

    positivity = fit.rows(output/'joint_positivity.csv')
    fig, axes = plt.subplots(2, 1, figsize=(10, 7))
    for s, name in enumerate(names):
        selected = [r for r in positivity if int(r['setting_index']) == s]
        x = [int(r['truth_block']) for r in selected]
        axes[0].plot(x, [float(r['minimum_response_bracket'])*scale for r in selected], 'o-', label=name)
        axes[1].plot(x, [float(r['epsilon_max']) for r in selected], 'o-', label=name)
        axes[1].axhline(nominal[s]['epsilon'], ls='--', color=f'C{s}', alpha=.6)
    axes[0].axhline(0, color='black', lw=.7)
    axes[0].set(ylabel=f'Minimum angular bracket [{unit}]', title='Positivity envelope by setting, including guard bins')
    axes[1].set(xlabel='Truth block', ylabel=r'$\epsilon_{max}$ (solid); nominal (dashed)')
    axes[0].legend(fontsize=9)
    save(fig, 'positivity_and_epsilon')
    return pages


def render(output, partons=False, warmups=10000, calls=100000):
    output = Path(output).resolve()
    manifest_path = output/'joint_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    fit.need(warmups > 0 and calls > 0, 'PARTONS warmups and calls must be positive')
    report = output/'plots'/'plot_manifest.json'
    report.parent.mkdir(exist_ok=True)
    if report.is_file():
        previous = json.loads(report.read_text())
        fit.need(previous.get('status') != 'complete',
                 f'plot report already complete: {output/"all_joint_xsec_plots.pdf"}')
        fit.need(previous.get('partons') == partons,
                 'retry an incomplete report with the same PARTONS option to avoid mixing plot sets')
    state = {'status': 'running', 'partons': partons, 'individual_fits_run': False,
             'residual_corrected_points': 'central_only_joint_influence_errors_not_computed',
             'boundary_errors': 'unavailable_no_individual_refit_toys',
             'partons_warmups': warmups, 'partons_calls': calls}
    report.write_text(json.dumps(state, indent=2)+'\n')
    try:
        pages = native_plots(output, manifest, partons, warmups, calls)
        pages += joint_figures(output, manifest)
        combined = output/'all_joint_xsec_plots.pdf'
        # pdfunite does not overwrite its input pages; combine individual pages
        # only, avoiding duplicate pages from per-setting report PDFs.
        subprocess.run(['pdfunite', *(str(p) for p in pages), str(combined)], check=True)
        state.update(status='complete', pages=[str(p.relative_to(output)) for p in pages],
                     page_count=len(pages), combined_pdf=combined.name)
        manifest['joint_partons_projection'] = 'per_setting_native_GK06_GPDGK19' if partons else 'not_requested'
        manifest['plots'] = {'status': 'complete', 'manifest': 'plots/plot_manifest.json', 'combined_pdf': combined.name}
        manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    except Exception as error:
        state.update(status='failed', error=str(error))
        raise
    finally:
        report.write_text(json.dumps(state, indent=2)+'\n')
    print(f"[PLOTS] {len(pages)} pages: {combined}", flush=True)
