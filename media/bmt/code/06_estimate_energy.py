#!/usr/bin/env python3
"""
Rough energy (kcal/day) estimate for a diet table:  mmol/day x molecular weight x
Atwater factor (carbohydrate/protein 4, fat 9, ethanol 7, butyrate 6, sorbitol 2.4,
mannitol 1.6 kcal/g).  Only compounds with a known MW and factor are counted
(micronutrients, minerals and non-starch polymeric fibres are ignored).
Starch exchanges are polymers; residue counts are read from the AGORA names:
strch1 = 11-12 glucose residues, strch2 = 3-4.

Usage
-----
    python 06_estimate_energy.py ../data/inputs/custom_diet_fluxes.csv  ../results/comparison_all_diets/all_media_long.csv

Part 1 treats the first CSV as a VMH Diet Designer export (mmol/human/day;
duplicates summed; D-amino-acid duplicates reported separately).
Part 2 prints the same estimate for every medium in all_media_long.csv AS STORED.
For the MICOM media this is NOT an intake: stored values were diluted (x0.2 or
x0.1) for compounds absorbed in the small intestine, so non-starch energy is
under-counted by roughly 1/dilution.  Use it only as a sanity check on scale.
"""
import sys
import pandas as pd
GLC=162.14; H2O=18.02
def starch(n): return n*GLC+H2O
MW={ # g/mol
 'glc_D':180.16,'fru':180.16,'gal':180.16,'sucr':342.30,'malt':342.30,'lcts':342.30,'cellb':342.30,'tre':342.30,'melib':342.3,'raffin':504.4,
 'ala_L':89.09,'ala_D':89.09,'arg_L':174.20,'asn_L':132.12,'asp_L':133.10,'asp_D':133.10,'cys_L':121.16,'gln_L':146.15,'glu_L':147.13,'gly':75.07,'his_L':155.16,
 'ile_L':131.17,'leu_L':131.17,'lys_L':146.19,'met_L':149.21,'phe_L':165.19,'pro_L':115.13,'pro_D':115.13,'ser_L':105.09,'thr_L':119.12,'trp_L':204.23,
 'tyr_L':181.19,'val_L':117.15,'Lcystin':240.30,'4hpro':131.13,'glyb':117.15,'orn':132.16,
 'ocdcea':282.47,'hdca':256.42,'ocdca':284.48,'lnlc':280.45,'lnlnca':278.43,'ttdca':228.37,'hdcea':254.41,'ddca':200.32,'arach':312.53,'octa':144.21,'hca':116.16,
 'dca':172.26,'ptdca':242.40,'hpdca':270.45,'vacc':282.47,'gndc':310.50,'arachd':304.5,'epa':302.45,'ercc':338.57,'nrvnc':366.6,'crvnc':328.5,'clpnd':330.5,'dcsa':340.6,
 'tecsa':368.6,'strdnc':276.4,'dlnlc':308.5,'dlnlncg':306.5,'lnlncg':278.4,'mrgrl':268.4,'hpdca10Z':268.4,'ttdcea':226.4,'ocdcea_':282.47,'chsterol':386.65,
 'but':88.11,'glyc':92.09,'etoh':46.07,'sbt_D':182.17,'mnl':182.17,'xylt':152.15,
}
KCAL={}  # kcal/g
for k in ['glc_D','fru','gal','sucr','malt','lcts','cellb','tre','melib','raffin']: KCAL[k]=('carb',4)
for k in ['ala_L','ala_D','arg_L','asn_L','asp_L','asp_D','cys_L','gln_L','glu_L','gly','his_L','ile_L','leu_L','lys_L','met_L','phe_L','pro_L','pro_D','ser_L','thr_L','trp_L','tyr_L','val_L','Lcystin','4hpro','glyb','orn']: KCAL[k]=('protein',4)
for k in ['ocdcea','hdca','ocdca','lnlc','lnlnca','ttdca','hdcea','ddca','arach','octa','hca','dca','ptdca','hpdca','vacc','gndc','arachd','epa','ercc','nrvnc','crvnc','clpnd','dcsa','tecsa','strdnc','dlnlc','dlnlncg','lnlncg','mrgrl','hpdca10Z','ttdcea','chsterol','glyc']: KCAL[k]=('fat',9)
KCAL['etoh']=('alcohol',7); KCAL['but']=('SCFA',6); KCAL['sbt_D']=('polyol',2.4); KCAL['mnl']=('polyol',1.6); KCAL['xylt']=('polyol',2.4)
def energy(flux_per_day, n1, n2):
    mw=dict(MW); mw['strch1']=starch(n1); mw['strch2']=starch(n2); k=dict(KCAL); k['strch1']=('starch',4); k['strch2']=('starch',4)
    rows=[]
    for m,v in flux_per_day.items():
        if v<=0 or m not in mw or m not in k: continue
        g=v*mw[m]/1000; rows.append((m,k[m][0],v,g,g*k[m][1]))
    return pd.DataFrame(rows,columns=['met','type','mmol_day','g_day','kcal_day'])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    d = pd.read_csv(sys.argv[1]); d.columns = ["r", "f"]
    d["m"] = d.r.str.replace(r"^EX_", "", regex=True).str.replace(r"\[e\]$", "", regex=True)
    al = {"oleic": "ocdcea", "palmtl": "hdcea", "tedca": "ttdca", "xol": "chsterol", "occa": "octa", "4hpro_LT": "4hpro"}
    d["m"] = d.m.replace(al); g = d.groupby("m").f.sum()
    for (n1, n2) in [(11, 3), (12, 4)]:
        e = energy(g, n1, n2); t = e.groupby("type")[["g_day", "kcal_day"]].sum()
        print(f"\n=== custom diet, strch1={n1}, strch2={n2} residues: {e.kcal_day.sum():.0f} kcal/day, {e.g_day.sum():.0f} g/day")
        print(t.round(0).sort_values("kcal_day", ascending=False).T.to_string())
        e2 = e[~e.met.isin(["ala_D", "asp_D", "pro_D"])]
        print("  excluding D-amino-acid duplicates:", round(e2.kcal_day.sum()), "kcal")
    L = pd.read_csv(sys.argv[2])
    print("\n=== every medium as stored (see docstring caveat); starch residues 12/4")
    for dd in L.diet.unique():
        x = L[L.diet == dd].set_index("metabolite").flux
        e = energy(x, 12, 4); t = e.groupby("type").kcal_day.sum(); st = t.get("starch", 0)
        print(f"{dd:30s} total {t.sum():8.0f}   starch {st:7.0f}   non-starch {t.sum() - st:8.0f}")
