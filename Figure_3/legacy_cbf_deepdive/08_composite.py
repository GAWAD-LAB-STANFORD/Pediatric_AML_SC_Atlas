#!/usr/bin/env python
"""Figure 3 composite (CBF deep-dive): A CBF states | B focus KM ; C targetable |
D HOX | E CBF GO. Rows fill width (aspect preserved, no crop). Letters only."""
import os, sys
from PIL import Image, ImageDraw, ImageFont
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
P = C.OUT
FILES = {"A":"panelA_cbf_states.png","B":"panelB_focus_km.png","C":"panelC_targetable.png",
         "D":"panelD_hox.png","E":"panelE_cbf_go.png"}
ROWS = [["A","B"], ["C","D","E"]]
W, GAP, HEAD, TOP, MARGIN = 2500, 46, 76, 26, 30; MAXH = 0.60*W
def font(sz):
    for p in ["/System/Library/Fonts/Supplemental/Arial Bold.ttf","/System/Library/Fonts/Helvetica.ttc"]:
        try: return ImageFont.truetype(p, sz)
        except Exception: pass
    return ImageFont.load_default()
fL = font(62)
missing = [k for k,v in FILES.items() if not os.path.exists(os.path.join(P,v))]
if missing: raise SystemExit("missing panels: "+", ".join(missing)+" (run steps 01-07 first)")
imgs = {k: Image.open(os.path.join(P, v)).convert("RGB") for k, v in FILES.items()}
rowlays=[]
for row in ROWS:
    ratios=[imgs[k].width/imgs[k].height for k in row]; n=len(row); avail=W-GAP*(n-1)
    Hrow=min(avail/sum(ratios), MAXH); widths=[r*Hrow for r in ratios]
    tot=sum(widths)+GAP*(n-1); x=MARGIN+(W-tot)/2; placed=[]
    for k,w in zip(row,widths): placed.append((k,int(round(x)),int(round(w)))); x+=w+GAP
    rowlays.append((int(round(Hrow)),placed))
totalH=TOP+sum(HEAD+h for h,_ in rowlays)+GAP*len(rowlays)
canvas=Image.new("RGB",(W+2*MARGIN,int(totalH)),"white"); d=ImageDraw.Draw(canvas); y=TOP
for Hrow,placed in rowlays:
    for k,x,w in placed:
        d.text((x,y),k,font=fL,fill="black"); canvas.paste(imgs[k].resize((w,Hrow),Image.LANCZOS),(x,y+HEAD))
    y+=HEAD+Hrow+GAP
out=os.path.join(os.path.dirname(P),"Figure_3__CBF_states")
canvas.save(out+".png"); canvas.save(out+".pdf","PDF",resolution=200)
print("wrote %s.png (%dx%d, aspect %.2f)"%(out,canvas.width,canvas.height,canvas.height/canvas.width))
