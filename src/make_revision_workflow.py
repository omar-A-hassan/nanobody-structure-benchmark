"""Draw the study workflow from verified cohort counts."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/revision_2026-10-23/publication'
OUT.mkdir(parents=True,exist_ok=True)
counts = json.loads((ROOT/'outputs/revision_2026-10-17/overlap_extension/overlap_extension_summary.json').read_text())
assert counts['cohorts']['original_common']['n']==588
assert counts['cohorts']['extended_lower_risk']['n']==393
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(9,6.2))
ax.set(xlim=(0,10),ylim=(0,8)); ax.axis('off')

def box(x,y,w,h,label,color='#eef3f6'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.06,rounding_size=0.05',
                              linewidth=1,edgecolor='#455a64',facecolor=color))
    ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=10.5,linespacing=1.4)

def arrow(start,end):
    ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>',mutation_scale=12,lw=1.1,color='#455a64'))

ax.text(0,7.8,'(a) Structural benchmark and SimpleFold analyses',weight='bold',fontsize=12)
box(.15,6.55,2.4,.8,'597 experimental VHHs\nSAbDab2 references')
box(3.2,6.55,6.4,.8,'VHH sequence only at inference\nSimpleFold 597   |   IgFold 597   |   NBB2 588')
arrow((2.65,6.95),(3.1,6.95))
box(.15,4.8,4.35,1.05,'Common structural comparison: 588\nH1, H2, H3 and framework\nWhole-chain and framework fits')
box(5.05,4.8,4.55,1.05,'Training-overlap and redundancy screens\n393 representatives\nAvailable-source sensitivity analysis')
arrow((4.9,6.45),(2.35,5.95));arrow((4.6,5.33),(4.95,5.33))
box(.15,3.1,3,.95,'Sampling disagreement\n597 targets\nLength-controlled error analysis')
box(3.55,3.1,2.8,.95,'Finite sampling coverage\n30 difficult + 10 controls\n20 samples × 4 settings')
box(6.75,3.1,2.85,.95,'ESM2 class probe\n518 labelled targets\nGrouped five-fold analysis')
# Route the independent SimpleFold analyses around, not through, the boxes.
ax.plot([9.68,9.9,9.9,1.65],[6.95,6.95,4.35,4.35],color='#455a64',lw=1.1)
for x in [1.65,4.95,8.15]: arrow((x,4.35),(x,4.15))
ax.plot([0,10],[2.65,2.65],color='#bcc7cd',lw=.8)
ax.text(0,2.3,'(b) Separate illustrative docking case: one complex',weight='bold',fontsize=12)
box(.15,.6,2.9,1.1,'Experimental VHH control\nand three predicted VHHs\nfrom the 7R73 case')
box(3.55,.6,2.65,1.1,'Fixed bound antigen\nSame HDOCK protocol\nTop 10 poses per input')
box(6.75,.6,2.85,1.1,'DockQ and components\nEpitope recovery\nDescriptive case results')
arrow((3.15,1.15),(3.45,1.15));arrow((6.3,1.15),(6.65,1.15))
fig.subplots_adjust(left=.02,right=.99,bottom=.02,top=.98)
for suffix in ['png','pdf','svg']:fig.savefig(OUT/f'workflow.{suffix}',dpi=240,bbox_inches='tight',pad_inches=.08)
print(OUT/'workflow.pdf')
