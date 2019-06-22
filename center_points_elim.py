


''' program starts here'''

'''resin rich area '''
wres = w-7*r
hres = h-7*r    
xres = random.random()*wres-wres/2
yres = random.random()*hres-hres/2
resc1 = [xres,yres]            
res_rad1 = random.random()*3*r+4*r#first one is between 2-4 times radiuses
#print res_rad1, resc1

xres_shift = res_rad1*random.random()*pi
yres_shift = res_rad1*random.random()*pi
resc2 = [xres+ xres_shift * 0.7 ,yres+ yres_shift*0.7]
res_rad2 = random.random()*3*r+4*r
#print res_rad2, resc2

removals = []
for cp_candid in all_cps:
    for cpc in cp_candid:
        d1 = ((cpc[0]-resc1[0])**2 + (cpc[1]-resc1[1])**2 ) ** 0.5
        d2 = ((cpc[0]-resc2[0])**2 + (cpc[1]-resc2[1])**2 ) ** 0.5
        if d1 < res_rad1 or d2 < res_rad2:
            removals.append(cp_candid)#        
            continue
        
        
for rem in removals:
    all_cps.remove(rem)
        









