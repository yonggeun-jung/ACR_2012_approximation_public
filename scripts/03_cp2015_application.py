"""
04_cp_application.py
ACR error bounds on the Caliendo and Parro (2015) NAFTA welfare calculation.

Inputs
------
From CP (2015) replication package
  (download from https://academic.oup.com/restud/article/82/1/1/1547758):
    ../data/cp2015/T.txt              Sectoral trade elasticities (20x1)
    ../data/cp2015/xbilat1993.txt     Bilateral trade flows 1993 (620x31)
    ../data/cp2015/alphas.csv         Final demand shares (40x31)
    ../data/cp2015/dlnpi_nafta.csv    Counterfactual Δln π^j_nn (40x31)

  The last two files are exported from CP's .mat workspace in MATLAB:
    load('initial_condition_1993_noS.mat'); load('alphas.mat');
    csvwrite('alphas.csv', alphas);
    % After running CP_counterfactuals.m (NAFTA-only block):
    for j=1:J
      pi_nn_base(j,:)=diag(Din(1+(j-1)*N:j*N,:));
      pi_nn_cf(j,:)=diag(Dinp_all_oN(1+(j-1)*N:j*N,:));
    end
    csvwrite('dlnpi_nafta.csv', log(pi_nn_cf)-log(pi_nn_base));

From CEPII Gravity (see scripts/02_build_dist.do; reads ../data/cepii/raw/Gravity_V202211.dta):
    ../data/cp2015/dist_cp31.csv

Outputs
-------
    ../output/tables/table_cp_gravity.csv
    ../output/tables/table_cp_bounds.csv
    ../output/tables/table_cp_bounds_summary.csv
"""

import numpy as np
import csv
import os
import sys
import statsmodels.api as sm
from statsmodels.genmod.families import Poisson

SCRIPT_DIR = os.path.dirname(__file__)
if SCRIPT_DIR:
    os.chdir(SCRIPT_DIR)

CP_DIR = os.path.join('..', 'data', 'cp2015')
OUT_DIR = os.path.join('..', 'output', 'tables')
os.makedirs(OUT_DIR, exist_ok=True)

# Constants
J, N = 20, 31
sigma = 4
MEX, CAN, US = 19, 4, 29

CP_COUNTRIES = [
    'ARG','AUS','AUT','BRA','CAN','CHL','CHN','DNK','FIN','FRA',
    'DEU','GRC','HUN','IND','IDN','IRL','ITA','JPN','KOR','MEX',
    'NLD','NZL','NOR','PRT','ZAF','ESP','SWE','TUR','GBR','USA','ROW'
]
CP_SECTORS = [
    'Agriculture','Mining','Food','Textile','Wood','Paper','Petroleum',
    'Chemicals','Plastic','Minerals','Basic metals','Metal products',
    'Machinery nec','Office','Electrical','Communication','Medical',
    'Auto','Other Transport','Other Manuf'
]

# Load data
T       = np.loadtxt(os.path.join(CP_DIR, 'T.txt'))
xbilat  = np.loadtxt(os.path.join(CP_DIR, 'xbilat1993.txt'))
alphas  = np.loadtxt(os.path.join(CP_DIR, 'alphas.csv'), delimiter=',')
dlnpi   = np.loadtxt(os.path.join(CP_DIR, 'dlnpi_nafta.csv'), delimiter=',')

# Distances (30 countries, excluding ROW)
iso3_to_idx = {c: i for i, c in enumerate(CP_COUNTRIES)}
dist_info = {}
with open(os.path.join(CP_DIR, 'dist_cp31.csv')) as f:
    for row in csv.DictReader(f):
        o, d = row['iso3_o'], row['iso3_d']
        if o in iso3_to_idx and d in iso3_to_idx:
            n, i = iso3_to_idx[d], iso3_to_idx[o]
            dist_info[(n, i)] = {
                'dist':   float(row['dist']),
                'contig': int(row['contig']),
                'comlang': int(row['comlang_off']),
                'colony': int(row['col_dep_ever']),
            }

# Build regression panel (common across sectors)
pairs = [(n, i) for n in range(30) for i in range(30)
         if n != i and (n, i) in dist_info]
K = len(pairs)

ln_d  = np.array([np.log(dist_info[p]['dist']) for p in pairs])
ln_d2 = ln_d ** 2
ctrl  = np.column_stack([
    np.array([dist_info[p]['contig'] for p in pairs]),
    np.array([dist_info[p]['comlang'] for p in pairs]),
    np.array([dist_info[p]['colony'] for p in pairs]),
])

exp_ids = sorted({p[1] for p in pairs})
imp_ids = sorted({p[0] for p in pairs})
exp_dum = np.zeros((K, len(exp_ids) - 1))
imp_dum = np.zeros((K, len(imp_ids) - 1))
for k, (n, i) in enumerate(pairs):
    ei = exp_ids.index(i)
    if ei > 0: exp_dum[k, ei - 1] = 1
    mi = imp_ids.index(n)
    if mi > 0: imp_dum[k, mi - 1] = 1

ln_d_min, ln_d_max = ln_d.min(), ln_d.max()

X_quad = sm.add_constant(np.column_stack([ln_d, ln_d2, ctrl, exp_dum, imp_dum]))

# Sectoral gravity diagnostic
diag = {}
for j in range(J):
    y = np.array([xbilat[j * N + n, i] for (n, i) in pairs])
    try:
        res = sm.GLM(y, X_quad, family=Poisson()).fit(disp=False, cov_type='HC1')
        b2    = res.params[2]
        se_b2 = res.bse[2]
        p_b2  = res.pvalues[2]
        deps  = 2 * abs(b2) * (ln_d_max - ln_d_min)
    except Exception:
        b2, se_b2, p_b2, deps = 0., 0., 1., 0.
    diag[j] = {'beta2': b2, 'se_b2': se_b2, 'p': p_b2, 'delta_eps': deps}

# Export gravity results
grav_path = os.path.join(OUT_DIR, 'table_cp_gravity.csv')
with open(grav_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['sector','theta','beta2','se_b2','p_b2','delta_eps'])
    for j in range(J):
        d = diag[j]
        w.writerow([CP_SECTORS[j], f"{T[j]:.2f}",
                    f"{d['beta2']:.6f}", f"{d['se_b2']:.6f}",
                    f"{d['p']:.6f}", f"{d['delta_eps']:.4f}"])

# Compute bounds
rows = []
totals = {c: 0. for c in ['Mexico','Canada','US']}
totals_w = {c: 0. for c in ['Mexico','Canada','US']}

for j in range(J):
    de = diag[j]['delta_eps']
    row = {'sector': CP_SECTORS[j], 'theta': T[j], 'delta_eps': de,
           'beta2': diag[j]['beta2'], 'p': diag[j]['p']}
    for c_idx, cn in [(MEX,'Mexico'),(CAN,'Canada'),(US,'US')]:
        v = dlnpi[j, c_idx]
        dp = abs(v) if not np.isnan(v) else 0.
        aj = alphas[j, c_idx]
        bnd  = (de / (sigma - 1)) * dp * 100
        wbnd = (de / (sigma - 1)) * dp * aj / T[j] * 100
        wc   = aj / T[j] * dp * 100
        row[f'dlnpi_{cn}'] = dp
        row[f'alpha_{cn}'] = aj
        row[f'bnd_{cn}'] = bnd
        row[f'wbnd_{cn}'] = wbnd
        row[f'wcontr_{cn}'] = wc
        totals[cn] += wbnd
        totals_w[cn] += wc
    rows.append(row)

# Export bounds
bnd_path = os.path.join(OUT_DIR, 'table_cp_bounds.csv')
with open(bnd_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['sector','theta','delta_eps','beta2','p',
                'dlnpi_Mexico','alpha_Mexico','bnd_Mexico','wbnd_Mexico','wcontr_Mexico',
                'dlnpi_Canada','alpha_Canada','bnd_Canada','wbnd_Canada','wcontr_Canada',
                'dlnpi_US','alpha_US','bnd_US','wbnd_US','wcontr_US'])
    for r in rows:
        w.writerow([r['sector'], f"{r['theta']:.2f}",
                    f"{r['delta_eps']:.4f}", f"{r['beta2']:.6f}", f"{r['p']:.6f}",
                    f"{r['dlnpi_Mexico']:.6f}", f"{r['alpha_Mexico']:.4f}",
                    f"{r['bnd_Mexico']:.4f}", f"{r['wbnd_Mexico']:.6f}", f"{r['wcontr_Mexico']:.6f}",
                    f"{r['dlnpi_Canada']:.6f}", f"{r['alpha_Canada']:.4f}",
                    f"{r['bnd_Canada']:.4f}", f"{r['wbnd_Canada']:.6f}", f"{r['wcontr_Canada']:.6f}",
                    f"{r['dlnpi_US']:.6f}", f"{r['alpha_US']:.4f}",
                    f"{r['bnd_US']:.4f}", f"{r['wbnd_US']:.6f}", f"{r['wcontr_US']:.6f}"])

cp_w = {'Mexico': 1.31, 'Canada': -0.06, 'US': 0.08}
sum_path = os.path.join(OUT_DIR, 'table_cp_bounds_summary.csv')
with open(sum_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['country','cp_welfare_pp','bound_ek_pp'])
    for cn in ['Mexico','Canada','US']:
        w.writerow([cn, f"{cp_w[cn]:.2f}", f"{totals[cn]:.4f}"])

# Print summary
n5  = sum(1 for j in range(J) if diag[j]['p'] < 0.05)
n10 = sum(1 for j in range(J) if diag[j]['p'] < 0.10)
print(f"Sectors significant at 5%: {n5}/{J}")
print(f"Sectors significant at 10%: {n10}/{J}")
for cn in ['Mexico','Canada','US']:
    print(f"  {cn}: CP welfare = {cp_w[cn]:.2f} pp, EK bound = {totals[cn]:.4f} pp")
print('[Done]')