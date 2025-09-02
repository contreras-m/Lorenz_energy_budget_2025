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
        'BKE': (1, 7),
        'SMKE': (8, 7),
        'SMPE': (8,3),
        'BPE': (1, 3),
    }

    for label, (x, y) in boxes.items():
        ax.add_patch(plt.Rectangle((x, y), dimx_box, dimy_box, fill=True, edgecolor='navy',facecolor='navy'))
        ax.text(x + 1.5, y + 0.5, label, ha='center', va='center', fontsize=szfont,color='white',fontweight='bold')

    # Define arrows and labels for conversions
    def draw_doublearrow(start, end,label,dimx_box,dimy_box, side="right", offset=(0,0), sizefont=16):
        sx, sy = boxes[start]
        ex, ey = boxes[end]
        if side=="right":
            sx += dimx_box  # arrow start at right edge
            ey += dimy_box/2 
            sy += dimy_box/2
            ax.annotate('', xy=((sx+ex)/2-0.1, ey), xytext=(sx, sy),
            arrowprops=dict(arrowstyle='-|>',facecolor='black', lw=2,mutation_scale=20))
            midx = (sx + (sx+ex)/2) / 2 + offset[0]
            midy = (sy + ey) / 2 + offset[1]
            ax.text(midx, midy, label,fontsize=sizefont, ha='center', va='center',
            bbox=dict(facecolor='lightblue', edgecolor='black', boxstyle='round,pad=0.3'))
        elif side=="left": 
            sx += dimx_box  # arrow start at right edge
            ey += dimy_box/2 
            sy += dimy_box/2
            ax.annotate('', xy=(ex, ey), xytext=((sx+ex)/2+0.1, sy),
            arrowprops=dict(arrowstyle='-|>',facecolor='black', lw=2,mutation_scale=20))
            midx = (sx + (sx+ex)) / 2 + offset[0]
            midy = (sy + ey) / 2 + offset[1]
            ax.text(midx, midy, label,fontsize=sizefont, ha='center', va='center',
            bbox=dict(facecolor='lightblue', edgecolor='black', boxstyle='round,pad=0.3'))
        
        
        
    # Define arrows and labels for conversions
    def draw_arrow(start, end,label,dimx_box,dimy_box, direction="h", offset=(0,0), sizefont=16):
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
                    arrowprops=dict(arrowstyle='-|>',facecolor='black', lw=2,mutation_scale=20))
        midx = (sx + ex) / 2 + offset[0]
        midy = (sy + ey) / 2 + offset[1]
        ax.text(midx, midy, label,fontsize=sizefont, ha='center', va='center',
        bbox=dict(facecolor='lightblue', edgecolor='black', boxstyle='round,pad=0.3'))

    # Define arrows and labels for conversions
    def draw_varrow(start, end,label,dimx_box,dimy_box,sign="negative", offset=(0,0),sizefont=16):
        sx, sy = boxes[start]
        ex, ey = boxes[end]

        if sign=="negative":
            sy+=dimy_box/2
            ey+=dimy_box/2   
            ex += dimx_box/2 
            sx += dimx_box/2
            ax.annotate('', xy=((sx+ex)/2,sy+2 ), xytext=((sx+ex)/2, ey),
                    arrowprops=dict(arrowstyle='-|>',facecolor='black',lw=2,mutation_scale=20))
        elif sign=="positive":
            sy-=3*dimy_box/2
            ey-=3*dimy_box/2   
            ex += dimx_box/2 
            sx += dimx_box/2
            ax.annotate('', xy=((sx+ex)/2,sy ), xytext=((sx+ex)/2, ey+2),
                    arrowprops=dict(arrowstyle='-|>',facecolor='black',lw=2,mutation_scale=20))
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
    def draw_residual(start, label1,label2a,label2b,label2c,label3a,label3b,dimx_box,dimy_box,
                      ad2a='-|>',ad2b='-|>',ad2c='-|>',ad3a='-|>',ad3b='-|>',lines=3,offset=(0,0),sizefont=16):
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
        sxmini, symini = boxesmini[label1]
        if lines==3:
            boxesmini2 = {label2a: (sxmini-1.2, symini+1),label2b: (sxmini+0.5, symini+1),label2c: (sxmini+2.2, symini+1),}
        elif lines==2:
            boxesmini2 = {label2a: (sxmini-.4, symini-1),label2b: (sxmini+1.3, symini-1),}

        width=1.4
        height=0.3
        for label, (x, y) in boxesmini2.items():
            rounded_box = FancyBboxPatch(
            (x, y),
            width, height,
            boxstyle="round,pad=0.1,rounding_size=0.2",
            linewidth=1,
            edgecolor='black',
            facecolor='khaki'
            )
            ax.add_patch(rounded_box)
            ax.text(x+0.7 , y+0.15, label, ha='center', va='center', fontsize=sizefont-3)

        if lines==2:
            draw_arrow2(boxesmini[label1],boxesmini2[label2b], r'',1,0.45, direction="v",location="bottom",arrowdir=ad2b, offset=(0, 0),offset_arrow=(0.5,-0.5))
            draw_arrow2(boxesmini[label1],boxesmini2[label2a], r'',1,0.45, direction="v",location="bottom",arrowdir=ad2a, offset=(0, 0),offset_arrow=(0.5,-0.5))
        if lines==3:
            draw_arrow2(boxesmini[label1],boxesmini2[label2b], r'',1,0.45, direction="v",arrowdir=ad2b, offset=(0, 0),offset_arrow=(0.5,0.1))
            draw_arrow2(boxesmini[label1],boxesmini2[label2a], r'',1,0.45, direction="v",arrowdir=ad2a, offset=(0, 0),offset_arrow=(0.5,0.1))
            draw_arrow2(boxesmini[label1],boxesmini2[label2c], r'',1,0.45, direction="v",arrowdir=ad2c, offset=(0, 0),offset_arrow=(0.5,0.1))
            #--------------------
            sxmini2, symini2 = boxesmini2[label2b]
            boxesmini3 = {label3a: (sxmini2-.895, symini2+1),label3b: (sxmini2+0.95, symini2+1),}

            width=1.4
            height=0.3
            for label, (x, y) in boxesmini3.items():
                rounded_box = FancyBboxPatch(
                (x, y),
                width, height,
                boxstyle="round,pad=0.1,rounding_size=0.2",
                linewidth=1,
                edgecolor='black',
                facecolor='beige'
                )
                ax.add_patch(rounded_box)
                ax.text(x+0.7 , y+0.15, label, ha='center', va='center', fontsize=sizefont-5)

            draw_arrow2(boxesmini2[label2b],boxesmini3[label3a], r'',1,0.4, direction="v",arrowdir=ad3a, offset=(0, 0),offset_arrow=(0,-0))
            draw_arrow2(boxesmini2[label2b],boxesmini3[label3b], r'',1,0.4, direction="v",arrowdir=ad3b, offset=(0, 0),offset_arrow=(0,-0))


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
    draw_doublearrow('BKE', 'SMKE', r'$CK_{(SM,B)}$' + '\n' + fr'$({mvar[3]:.2f})$',dimx_box,dimy_box, side="right", offset=(-0.1, 0),sizefont=12)
    draw_doublearrow('BKE', 'SMKE', r'$CK_{(B,SM)}$' + '\n' + fr'$( {mvar[2]:.2f})$',dimx_box,dimy_box, side="left", offset=(-1, 0),sizefont=12)
    draw_varrow('BKE', 'SMKE', r'$BnK_{(B,SM)}$' + '\n' + fr'$({mvar[3]+mvar[2]:.2f})$',dimx_box,dimy_box, offset=(0, 1),sizefont=12)

    draw_doublearrow('BPE', 'SMPE', r'$CP_{(SM,B)}$' + '\n' + fr'$({mvar[10]:.2f})$',dimx_box,dimy_box, side="right", offset=(-0.1, 0),sizefont=12)
    draw_doublearrow('BPE', 'SMPE', r'$CP_{(B,SM)}$' + '\n' + fr'$( {mvar[9]:.2f})$',dimx_box,dimy_box, side="left", offset=(-1, 0),sizefont=12)
    draw_varrow('BPE', 'SMPE', r'$BnP_{(B,SM)}$' + '\n' + fr'$({+mvar[10]+mvar[9]:.2f})$',dimx_box,dimy_box,sign="positive", offset=(0, 1),sizefont=12)

    draw_arrow('SMPE', 'SMKE', fr'$PK_{{SM}} \ ({mvar[6]:.2f})$' ,dimx_box,dimy_box, direction="v", offset=(0, 0),sizefont=18)
    draw_arrow('BPE', 'BKE', fr'$PK_{{B}} \ ({mvar[11]:.2f})$',dimx_box,dimy_box, direction="v", offset=(0, 0),sizefont=18)

    draw_residual('BKE',r'$FK_B+DK_B$' + fr'$({mvar[24]:.2f})$',r'$\tau(b)_B$' +fr'$({-mvar[14]:.2f})$' ,r'$\tau(s)_B$' +fr'$({mvar[16]:.2f})$',
                  r'$DK_B$'+ fr'$({mvar[28]:.2f})$',r'$\tau(s)_B^{geos}$'+fr'$({mvar[20]:.2f})$',r'$\tau(s)_B^{ageos}$'+fr'$({mvar[18]:.2f})$',dimx_box,dimy_box,
                   ad2a='-|>',ad2b='<|-',ad2c='-|>',ad3a='<|-',ad3b='<|-',lines=3, offset=(0,0),sizefont=15)
    draw_residual('SMKE',r'$FK_{SM}+DK_{SM}$'+ fr'$({mvar[22]:.2f})$',r'$\tau(b)_{SM}$'+fr'$({-mvar[15]:.2f})$' ,r'$\tau(s)_{SM}$'+fr'$({mvar[17]:.2f})$' ,
                  r'$DK_{SM}$'+ fr'$({mvar[26]:.2f})$',r'$\tau(s)_{SM}^{geos}$'+fr'$({mvar[21]:.2f})$' ,r'$\tau(s)_{SM}^{ageos}$'+fr'$({mvar[19]:.2f})$' ,dimx_box,dimy_box,
                   ad2a='-|>',ad2b='<|-',ad2c='-|>',ad3a='-|>',ad3b='<|-',lines=3, offset=(0,0),sizefont=15)



    draw_onebox('SMKE', fr'$P_{{SM}} \ ({mvar[7]-mvar[6]:.2f})$',dimx_box,dimy_box,location="bottom",arrowdir='<-'  ,position=(-0.2, -.8), offset=(0.7, 0.15),sizefont=12)
    draw_onebox('SMKE', fr'$AK_{{SM}} \  ({mvar[0]:.2f})$',dimx_box,dimy_box,location="bottom",arrowdir='<-' ,position=(+1.8, -0.8), offset=(0.7, 0.15),sizefont=12)

    draw_onebox('BKE',fr'$P_{{B}} \ ({mvar[8]-mvar[11]:.2f})$',dimx_box,dimy_box,location="bottom" ,position=(-0.2, -.8), offset=(0.7, 0.15),sizefont=12)
    draw_onebox('BKE',fr'$AK_{{B}} \ ({mvar[1]:.2f})$',dimx_box,dimy_box,location="bottom",arrowdir='<-'  ,position=(+1.8, -0.8), offset=(0.7, 0.15),sizefont=12)

    draw_onebox('SMPE', fr'$AP_{{SM}} \  ({mvar[4]:.2f})$',dimx_box,dimy_box,location="upper",arrowdir='->'  ,position=(+1.8, 1.5), offset=(0.7, 0.15),sizefont=12)
    draw_onebox('BPE', fr'$AP_{{B}} \ ({mvar[5]:.2f})$',dimx_box,dimy_box,location="upper" ,arrowdir='<-',position=(+1.8, 1.5), offset=(0.7, 0.15),sizefont=12)

    draw_residual('BPE',fr'$FP_B+DP_B \ ({mvar[25]:.2f})$',r'$FP_B$'+ fr'$({mvar[12]:.2f})$',r'$DP_B$'+ fr'$({mvar[29]:.2f})$','','','',
                  dimx_box,dimy_box, ad2a='<|-',ad2b='-|>',lines=2, offset=(0,0),sizefont=16)
    draw_residual('SMPE',fr'$FP_{{SM}}+DP_{{SM}} \ ({mvar[23]:.2f})$',r'$FP_{{SM}}$'+ fr'$({mvar[13]:.2f})$',r'$DP_{{SM}}$'+ fr'$({mvar[27]:.2f})$','','','',
                  dimx_box,dimy_box, ad2a='<|-',ad2b='<|-',lines=2, offset=(0,0),sizefont=16)

    plt.savefig('/Users/contrema/Documents/lorenz2.png',dpi=200, bbox_inches='tight')


    plt.show()


#plot_lorenz_energy_cycle(m)

