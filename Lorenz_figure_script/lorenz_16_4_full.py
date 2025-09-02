import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

def plot_lorenz_energy_cycle(mvar):
    dimx_box=3
    dimy_box=1
    szfont=18
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 12)
    ax.axis('off')

    # Define box positions and labels
    boxes = {
        r'$BKE^{16km}$': (1, 7),
        r'$SMKE^{16km}$': (8, 7),
        r'$SMPE^{16km}$': (8,3),
        r'$BPE^{16km}$': (1, 3),
    }

    for label, (x, y) in boxes.items():
        ax.add_patch(plt.Rectangle((x, y), dimx_box, dimy_box, fill=True, edgecolor='navy',facecolor='navy'))
        ax.text(x + 1.5, y + 0.5, label, ha='center', va='center', fontsize=szfont,color='white',fontweight='bold')

    # Define arrows and labels for conversions
    def draw_doublearrow(start, end,label,dimx_box,dimy_box,arrowdir='->', side="right", offset=(0,0), sizefont=16):
        sx, sy = boxes[start]
        ex, ey = boxes[end]
        if side=="right":
            sx += dimx_box  # arrow start at right edge
            ey += dimy_box/2 
            sy += dimy_box/2
            ax.annotate('', xy=((sx+ex)/2-0.1, ey), xytext=(sx, sy),
            arrowprops=dict(arrowstyle=arrowdir,facecolor='black', lw=2,mutation_scale=20))
            midx = (sx + (sx+ex)/2) / 2 + offset[0]
            midy = (sy + ey) / 2 + offset[1]
            ax.text(midx, midy, label,fontsize=sizefont, ha='center', va='center',
            bbox=dict(facecolor='lightblue', edgecolor='black', boxstyle='round,pad=0.3'))
        elif side=="left": 
            sx += dimx_box  # arrow start at right edge
            ey += dimy_box/2 
            sy += dimy_box/2
            ax.annotate('', xy=(ex, ey), xytext=((sx+ex)/2+0.1, sy),
            arrowprops=dict(arrowstyle=arrowdir,facecolor='black', lw=2,mutation_scale=20))
            midx = (sx + (sx+ex)) / 2 + offset[0]
            midy = (sy + ey) / 2 + offset[1]
            ax.text(midx, midy, label,fontsize=sizefont, ha='center', va='center',
            bbox=dict(facecolor='lightblue', edgecolor='black', boxstyle='round,pad=0.3'))
        
        
        
    # Define arrows and labels for conversions
    def draw_arrow(start, end,label,dimx_box,dimy_box, arrowdir='->', direction="h", offset=(0,0), sizefont=16):
        sx, sy = boxes[start]
        ex, ey = boxes[end]
        if direction=="h":
            sx += dimx_box  # arrow start at right edge
            ey += dimy_box/2 
            sy += dimy_box/2 
        elif direction=="v": 
            sy+=dimy_box            
            ex += dimx_box/2 
            sx += dimx_box/2 
        ax.annotate('', xy=(ex, ey), xytext=(sx, sy),
                    arrowprops=dict(arrowstyle=arrowdir,facecolor='black', lw=2,mutation_scale=20))
        midx = (sx + ex) / 2 + offset[0]
        midy = (sy + ey) / 2 + offset[1]
        ax.text(midx, midy, label,fontsize=sizefont, ha='center', va='center',
        bbox=dict(facecolor='lightblue', edgecolor='black', boxstyle='round,pad=0.3'))

    # Define arrows and labels for conversions
    def draw_varrow(start, end,label,dimx_box,dimy_box,arrowdir='->',sign="negative", offset=(0,0),sizefont=16):
        sx, sy = boxes[start]
        ex, ey = boxes[end]

        if sign=="negative":
            sy+=dimy_box/2
            ey+=dimy_box/2   
            ex += dimx_box/2 
            sx += dimx_box/2
            ax.annotate('', xy=((sx+ex)/2,sy+2 ), xytext=((sx+ex)/2, ey),
                    arrowprops=dict(arrowstyle=arrowdir,facecolor='black',lw=2,mutation_scale=20))
        elif sign=="positive":
            sy-=3*dimy_box/2
            ey-=3*dimy_box/2   
            ex += dimx_box/2 
            sx += dimx_box/2
            ax.annotate('', xy=((sx+ex)/2,sy ), xytext=((sx+ex)/2, ey+2),
                    arrowprops=dict(arrowstyle=arrowdir,facecolor='black',lw=2,mutation_scale=20))
        midx = (sx + ex) / 2 + offset[0]
        midy = (sy + ey) / 2 + offset[1]
        ax.text(midx, midy, label,fontsize=sizefont, ha='center', va='center',
        bbox=dict(facecolor='lightblue', edgecolor='black', boxstyle='round,pad=0.3'))


    # Define arrows and labels for conversions
    def draw_arrow2(start, end,label,dimx_box,dimy_box, direction="h",location="upper",arrowdir='->' , offset=(0,0),offset_arrow=(0,0)):
        sx, sy =start
        ex, ey = end


        if location=="upper":
            if direction=="h":
                sx += dimx_box  # arrow start at right edge
                ey += dimy_box/2 
                sy += dimy_box/2 
            elif direction=="v": 
                sy+=dimy_box           
                ex += dimx_box/2 
                sx += dimx_box/2
                ey -= dimy_box/4
            
            ax.annotate('', xy=(ex, ey), xytext=(sx+offset_arrow[0], sy+offset_arrow[1]),
            arrowprops=dict(arrowstyle=arrowdir,facecolor='black', lw=2))
        elif location=="bottom":
            if direction=="h":
                sx += dimx_box  # arrow start at right edge
                ey += dimy_box/2 
                sy += dimy_box/2 
            elif direction=="v": 
                sy+=dimy_box           
                ex += dimx_box/2 
                sx += dimx_box/2
                ey += dimy_box
            ax.annotate('', xy=(ex, ey), xytext=(sx+offset_arrow[0], sy+offset_arrow[1]),
            arrowprops=dict(arrowstyle=arrowdir,facecolor='black', lw=2))
       

       
        
    # Define arrows and labels for conversions
    def draw_residual(start, label1,dimx_box,dimy_box,
                      arrowdir='->',lines=3,offset=(0,0),sizefont=16):
        sx, sy = boxes[start]
        if lines==3:
            boxesmini = {label1: (sx+dimx_box/10, sy+dimx_box/2),}
            sx+=dimx_box/2
            sy+=dimy_box
            ax.annotate('', xy=(sx,sy), xytext=(sx, sy+.45),
            arrowprops=dict(arrowstyle='-|>',facecolor='black', lw=2))


        elif lines==2:
            boxesmini = {label1: (sx+dimx_box/10, sy-dimx_box/3),}
            sx+=dimx_box/2
            sy-=dimy_box
            ax.annotate('', xy=(sx,sy+1), xytext=(sx, sy+0.5),
            arrowprops=dict(arrowstyle='-|>',facecolor='black', lw=2))
      


      
        width=2.5
        height=0.45
        for label, (x, y) in boxesmini.items():
            rounded_box = FancyBboxPatch(
            (x, y),
            width, height,
            boxstyle="round,pad=0.1,rounding_size=0.2",
            linewidth=1,
            edgecolor='black',
            facecolor='goldenrod'
            )
            ax.add_patch(rounded_box)
            #ax.add_patch(plt.Rectangle((x, y), 2,0.6, fill=True, edgecolor='darkgray',facecolor='gray'))
            ax.text(x + 1.25, y+0.25 , label, ha='center', va='center', fontsize=sizefont)
        #--------------------
     
    def draw_onebox(start, label,dimx_box,dimy_box,position=(0,0), location="upper",arrowdir='->', offset=(0,0),sizefont=16):
        sx, sy = boxes[start]

        boxesmini = {label: (sx+position[0], sy+position[1]),}

        width=1.4
        height=0.3

        for label, (x, y) in boxesmini.items():
          rounded_box = FancyBboxPatch(
            (x, y),
            width, height,
            boxstyle="round,pad=0.1,rounding_size=0.2",
            linewidth=1,
            edgecolor='black',
            facecolor='aliceblue'
         )
          ax.add_patch(rounded_box)
          ax.text(x + offset[0], y+offset[1] , label, ha='center', va='center', fontsize=sizefont)
        sx, sy = boxesmini[label]
        sx+=dimx_box/4+0
        sy+=dimy_box/2
        
        if location=="bottom":
            ax.annotate('', xy=(sx,sy+0.35), xytext=(sx, sy-.1),
            arrowprops=dict(arrowstyle=arrowdir, lw=2))
        elif location=="upper":
            ax.annotate('', xy=(sx,sy-0.5), xytext=(sx, sy-1),
            arrowprops=dict(arrowstyle=arrowdir, lw=2))


        #midx = sx+ offset[0]
        #midy = sy  + offset[1]
       
        
    # Conversion arrows
    draw_doublearrow( r'$BKE^{16km}$',  r'$SMKE^{16km}$', r'$CK_{(SM,B)}^{16km}$' + '\n' + fr'$({mvar[3]:.2f})$',dimx_box,dimy_box,arrowdir='<-', side="right", offset=(-0.1, 0),sizefont=12)
    draw_doublearrow(r'$BKE^{16km}$', r'$SMKE^{16km}$', r'$CK_{(B,SM)}^{16km}$' + '\n' + fr'$( {mvar[2]:.2f})$',dimx_box,dimy_box,arrowdir='<-', side="left", offset=(-1, 0),sizefont=12)
    draw_varrow( r'$BKE^{16km}$',r'$SMKE^{16km}$', r'$BnK_{(B,SM)}^{16km}$' + '\n' + fr'$({mvar[3]+mvar[2]:.2f})$',dimx_box,dimy_box,arrowdir='->', offset=(0, 1),sizefont=12)

    draw_doublearrow(r'$BPE^{16km}$', r'$SMPE^{16km}$', r'$CP_{(SM,B)}^{16km}$' + '\n' + fr'$({mvar[10]:.2f})$',dimx_box,dimy_box, side="right", offset=(-0.1, 0),sizefont=12)
    draw_doublearrow(r'$BPE^{16km}$', r'$SMPE^{16km}$', r'$CP_{(B,SM)}^{16km}$' + '\n' + fr'$( {mvar[9]:.2f})$',dimx_box,dimy_box, side="left", offset=(-1, 0),sizefont=12)
    draw_varrow(r'$BPE^{16km}$', r'$SMPE^{16km}$', r'$BnP_{(B,SM)}^{16km}$' + '\n' + fr'$({+mvar[10]+mvar[9]:.2f})$',dimx_box,dimy_box,sign="positive", offset=(0, 1),sizefont=12)

    draw_arrow(r'$SMPE^{16km}$', r'$SMKE^{16km}$', fr'$PK_{{SM}}^{{16km}} \ ({mvar[6]:.2f})$' ,dimx_box,dimy_box,arrowdir='->', direction="v", offset=(0, 0),sizefont=18)
    draw_arrow( r'$BPE^{16km}$',  r'$BKE^{16km}$', fr'$PK_{{B}}^{{16km}} \ ({mvar[11]:.2f})$',dimx_box,dimy_box, direction="v",arrowdir='->', offset=(0, 0),sizefont=18)

    draw_residual( r'$BKE^{16km}$',r'$FK_B^{16km}+DK_B^{16km}$' + fr'$({mvar[24]:.2f})$',dimx_box,dimy_box,lines=3, offset=(0,0),sizefont=12)
    draw_residual(r'$SMKE^{16km}$',r'$FK_{SM}^{16km}+DK_{SM}^{16km}$'+ fr'$({mvar[22]:.2f})$',dimx_box,dimy_box,lines=3,arrowdir='->', offset=(0,0),sizefont=12)



    draw_onebox(r'$SMKE^{16km}$', fr'$P_{{SM}}^{{16km}} \ ({mvar[7]-mvar[6]:.2f})$',dimx_box,dimy_box,location="bottom",arrowdir='<-'  ,position=(-0.2, -.8), offset=(0.7, 0.15),sizefont=12)
    draw_onebox(r'$SMKE^{16km}$', fr'$AK_{{SM}}^{{16km}} \  ({mvar[0]:.2f})$',dimx_box,dimy_box,location="bottom",arrowdir='->' ,position=(+1.8, -0.8), offset=(0.7, 0.15),sizefont=12)

    draw_onebox(r'$BKE^{16km}$',fr'$P_{{B}}^{{16km}} \ ({mvar[8]-mvar[11]:.2f})$',dimx_box,dimy_box,location="bottom",arrowdir='<-'  ,position=(-0.2, -.8), offset=(0.7, 0.15),sizefont=12)
    draw_onebox(r'$BKE^{16km}$',fr'$AK_{{B}}^{{16km}} \ ({mvar[1]:.2f})$',dimx_box,dimy_box,location="bottom",arrowdir='<-'  ,position=(+1.8, -0.8), offset=(0.7, 0.15),sizefont=12)

    draw_onebox(r'$SMPE^{16km}$', fr'$AP_{{SM}}^{{16km}} \  ({mvar[4]:.2f})$',dimx_box,dimy_box,location="upper",arrowdir='<-'  ,position=(+1.8, 1.5), offset=(0.7, 0.15),sizefont=12)
    draw_onebox( r'$BPE^{16km}$', fr'$AP_{{B}}^{{16km}} \ ({mvar[5]:.2f})$',dimx_box,dimy_box,location="upper" ,arrowdir='<-',position=(+1.8, 1.5), offset=(0.7, 0.15),sizefont=12)

    draw_residual(r'$BPE^{16km}$',r'$FP_B^{16km}+DP_B^{16km}$' + fr'$({mvar[25]:.2f})$',dimx_box,dimy_box,lines=2, offset=(0,0),sizefont=12)
    draw_residual(r'$SMPE^{16km}$',r'$FP_{{SM}}^{16km}+DP_{{SM}}^{16km}$' + fr'$({mvar[23]:.2f})$',dimx_box,dimy_box,lines=2, offset=(0,0),sizefont=12)
    ax.text(1, 10 , '(a) 16 KM', ha='center', va='center', fontsize=22,fontweight='bold')

    
    plt.savefig('/Users/contrema/Documents/lorenz_16_full.png',dpi=200, bbox_inches='tight')
    

    plt.show()


#plot_lorenz_energy_cycle(m)

