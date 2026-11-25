# STRING figure (network without histones and hemoglobins): largest connected component,
# colored by direction (up/down), with the top 15 genes by degree. Output: A04_STRING_noHistHb.png/.pdf
import csv, networkx as nx, matplotlib, collections
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"DejaVu Sans"})
# Run from A_mass_trauma/04_network/string/  (requires networkx, matplotlib, adjustText)
B="../"
up={l.strip() for l in open(B+"input/A_DEG_up.txt") if l.strip()}
dn={l.strip() for l in open(B+"input/A_DEG_down.txt") if l.strip()}
ann=list(csv.reader(open(B+"string/DEG_reduced/string_protein_annotations_392_nodes.tsv"),delimiter="\t"))[1:]
dr={}
for r in ann:
    n={r[0]}|set(r[-1].split(",")); dr[r[0]]="up" if n&up else "down" if n&dn else "?"
E=list(csv.reader(open(B+"string/DEG_reduced/string_interactions_short_392_nodes.tsv"),delimiter="\t"))[1:]
G=nx.Graph(); G.add_edges_from((r[0],r[1]) for r in E)
H=G.subgraph(max(nx.connected_components(G),key=len)).copy()
deg=dict(G.degree())
col={"up":"#C0392B","down":"#2E6DA4","?":"#9E9E9E"}
fig,(ax,bx)=plt.subplots(1,2,figsize=(8.2,4.4),gridspec_kw={"width_ratios":[1.55,1]})
pos=nx.kamada_kawai_layout(H)
nx.draw_networkx_edges(H,pos,ax=ax,width=0.25,edge_color="#C8C8C8",alpha=0.55)
nodes=sorted(H,key=lambda g:deg[g])
nx.draw_networkx_nodes(H,pos,nodelist=nodes,ax=ax,node_size=[4+1.5*deg[g] for g in nodes],
    node_color=[col[dr[g]] for g in nodes],linewidths=0.3,edgecolors="white")
xs=[p[0] for p in pos.values()]; ys=[p[1] for p in pos.values()]
for g,(dx,dy) in {"KIF14":(-0.55,-0.55),"UBA52":(0.55,-0.45),"BCL2":(0.5,0.5),"COX5B":(-0.6,0.45)}.items():
    x,y=pos[g]
    ax.annotate(g,(x,y),xytext=(x+dx,y+dy),fontsize=8,fontweight="bold",color=col[dr[g]],
        ha="center",arrowprops=dict(arrowstyle="-",color="#424242",lw=0.6,shrinkA=2,shrinkB=4),
        bbox=dict(boxstyle="round,pad=0.2",fc="white",ec=col[dr[g]],lw=0.6))
ax.legend(handles=[Line2D([0],[0],marker="o",color="w",markerfacecolor=col["up"],markersize=7,label="Up-regulated in PTSD"),
    Line2D([0],[0],marker="o",color="w",markerfacecolor=col["down"],markersize=7,label="Down-regulated in PTSD")],
    loc="upper center",bbox_to_anchor=(0.5,0.02),ncol=2,frameon=False,fontsize=8)
ax.set_title("(a) Largest connected component",fontsize=9,loc="left")
ax.axis("off")
top=sorted(deg,key=lambda g:-deg[g])[:15][::-1]
bx.barh(top,[deg[g] for g in top],color=[col[dr[g]] for g in top],height=0.7)
for i,g in enumerate(top): bx.text(deg[g]+0.6,i,str(deg[g]),va="center",fontsize=7.5,color="#424242")
bx.set_xlabel("Node degree (STRING interactions)",fontsize=8.5)
bx.set_title("(b) Top 15 genes by degree",fontsize=9,loc="left")
bx.spines[["top","right"]].set_visible(False); bx.tick_params(labelsize=8)
bx.set_xlim(0,max(deg.values())*1.12)
plt.tight_layout(); plt.savefig("A04_STRING_noHistHb.png",dpi=300); plt.savefig("A04_STRING_noHistHb.pdf")
