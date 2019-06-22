def outskirt_check(cp_candid):
    xr = cp_candid[0]
    yr = cp_candid[1]
    out_center = []
    if xr < -w/2+r:
        xr1= w + xr
        out_center.append((xr1,yr))#, point1=(xr1+r, yr))
    if xr > w/2-r:
        xr1= -w + xr
        out_center.append((xr1,yr))
    if yr < -h/2 + r:
        yr1 = h + yr
        out_center.append((xr,yr1))
    if yr > h/2 - r:
        yr1 = -h + yr
        out_center.append((xr,yr1))
    # if the center is in the four corners
        
    if ( (xr + w/2)**2  + (yr +h/2)**2)**0.5 < r:
        xr1= w + xr
        yr1 = h + yr
        out_center.append((xr1,yr1))
    if ( (xr + w/2)**2 + (yr - h/2)**2)**0.5 < r:
        xr1= w + xr
        yr1 = -h + yr
        out_center.append((xr1,yr1))
    if ( (xr - w/2)**2 + (yr + h/2)**2)**0.5 < r:
        xr1= -w + xr
        yr1 = h + yr
        out_center.append((xr1,yr1))
    if ( (xr - w/2)**2 + (yr - h/2)**2)**0.5 < r:
        xr1= -w + xr
        yr1 = -h + yr
        out_center.append((xr1,yr1))

    return out_center

def distances(cp_candid, all_center_points):


    dist_cent = []
    for cpc in cp_candid:
        d1 = ((cpc[0]-resc1[0])**2 + (cpc[1]-resc1[1])**2 ) ** 0.5
        d2 = ((cpc[0]-resc2[0])**2 + (cpc[1]-resc2[1])**2 ) ** 0.5
        if d1 < res_rad1 or d2 < res_rad2:
            dist_cent.append(2*r-0.1)#putting a small value if the distance is smaller than resin radius
                
        for combined_point in all_center_points:
            for sp in combined_point:
                d = ((cpc[0]-sp[0])**2 + (cpc[1]-sp[1])**2 ) ** 0.5
                dist_cent.append(d)
           
    return dist_cent


def closest_cents(all_cps):
    cp = len(all_cps)
    main_points = []
    for i in range(len(all_cps)):
        main_points.append(all_cps[i][0])

    all_dist_cent = []
    for i in range(cp):
        copy_all_cps = copy.deepcopy(all_cps)#[:] #to prevent shallow copy!!!! (very BS feature)
        copy_all_cps[i] = [(h*h,w*w)] # removing that center from the rest and putting a large number
        mp = main_points[i]
        dist_cent_one = []
        for i in range(cp):
            dist_cent_one.append(min(distances([mp], [copy_all_cps[i]])))     # from all the points for one points choose the closest one for distance

        all_dist_cent.append(dist_cent_one) # distances from

    #finding n closest centers
    alldist = []
    alldist = copy.deepcopy(all_dist_cent)
    closest_centers = 4
    all_close_centers = []
    for i in range(cp):
        close_centers = []
        for j in range(closest_centers):
            ind_min = np.argmin(alldist[i])
            close_centers.append(ind_min)
            alldist[i][ind_min] = h*w #to replace the minimum
        all_close_centers.append(close_centers)

    sum_close_dist = [0]* cp#len(center_points)
    for i in range(cp):
        for j in all_close_centers[i]:
            sum_close_dist[i] += all_dist_cent[i][j]

    return sum_close_dist, all_close_centers
# this last part was finding the most isolated one
# the index of the maximum of that values returns the most isolated that is a candidate to move


def random_points(all_cps):
    cp = len(all_cps)
    mc = 0
    max_cycle = 10000
    while  mc < max_cycle: #and cp < fiber_num
        xr = random.random()*w-w/2
        yr = random.random()*h-h/2
        cp_candid = (xr,yr)
        center_points = [cp_candid]
        mc += 1
        out_center_candid = outskirt_check(cp_candid)
        if len(out_center_candid) > 0:
            for oc in out_center_candid:
                center_points.append(oc)

        dist_cent = distances(center_points, all_cps)
        if min(dist_cent) > 2*r+del_min:
            cp += 1
            mc = 0
            all_cps.append(center_points)
            
        if cp == max_cp:
            break

    return all_cps


def moving(first, second):
    move_vec = [(fi-se) for fi, se in zip(first,second)] #moving vector
    move_vec_len = (move_vec[0]**2  + move_vec[1]**2)**0.5
    ratio =  (move_vec_len - 2*r) / move_vec_len
    ratio = ratio* random.random()
    #del_min = 0.1*r
    if ratio < del_min:
        ratio =  del_min
        
    move_vec = [mv*ratio for mv in move_vec]

    new_point = [(fi-mv) for fi, mv in zip(first, move_vec)]  #adding vector to the moving candidate
    nps = [(new_point[0], new_point[1])]
    out_center_np = outskirt_check(new_point)
    if len(out_center_np) > 0:
        for oc in out_center_np:
            nps.append(oc)

    return nps




''' program starts here'''

'''resin rich area '''
wres = w-7*r
hres = h-7*r    
xres = random.random()*wres-wres/2
yres = random.random()*hres-hres/2
resc1 = [xres,yres]            
res_rad1 = 0.0#random.random()*3*r+4*r#first one is between 2-4 times radiuses
print res_rad1, resc1
xres_shift = res_rad1*random.random()*pi
yres_shift = res_rad1*random.random()*pi
resc2 = [xres+ xres_shift * 0.7 ,yres+ yres_shift*0.7]
res_rad2 = 0.0#random.random()*3*r+4*r
print res_rad2, resc2




xr = random.random()*w-w/2
yr = random.random()*h-h/2
cp_candid = (xr,yr) # adding first center point
center_points = [cp_candid]
out_center_candid = outskirt_check(cp_candid)
if len(out_center_candid) > 0:
    for oc in out_center_candid:
        center_points.append(oc)

all_cps = []
all_cps.append(center_points)

cp = 1
all_cps = random_points(all_cps)
cp = len(all_cps)

##center_points = [(0.0,0.0),(-10.0,-10.0),(10.0,-10.0),(10.0,10.0),(-10.0,10.0),
##(-5.0,-5.0),(5.0,5.0), (5.0,-5.0),(-5.0,5.0),
##(-0.0,-10.0), (0.0,10.0), (10.0,0.0),(-10.0,0.0) ]
##cp =  max(number_cp)
#print cp, vf


''' stirring the centerpoints around '''
# make moving a function moving the fiber, and the fiber that is moving too
max_add_steps = 120
add_step = 0
while add_step < max_add_steps:
    number_move = cp/2
    fiber_moves = 0
    trial_number = 12000#00
    trial =  0
    sum_close_dist, all_close_centers = closest_cents(all_cps)
    sum_close_dist_2 = copy.deepcopy(sum_close_dist)
    num_mov_can = cp/3
    mo_ca = []
    for i in range(num_mov_can): #number of moving candidate
        mci = np.argmax(sum_close_dist_2)
        mo_ca.append(mci)
        sum_close_dist_2[mci] = 0.0

    print mo_ca

    while fiber_moves < number_move and trial < trial_number:
        trial +=1
        k = trial % num_mov_can
        j = trial % len(all_close_centers[0])
        #print k
        #mo_ca = np.argmax(sum_close_dist) # moving candidate
        copy_all_cps = copy.deepcopy(all_cps)


        first = copy_all_cps[mo_ca[k]][0] #only main point
        second = copy_all_cps[all_close_centers[mo_ca[k]][j]][0] # changes from closest center to the others
        #move_vec = []
        nps = moving(first, second)



        rest_of_points = copy_all_cps[:mo_ca[k]] + copy_all_cps[mo_ca[k]+1:]
        d_np = distances(nps, rest_of_points)
        #del_min = 0.1*r
        if min(d_np) > 2*r+del_min:
            copy_all_cps = copy_all_cps[:mo_ca[k]] + [nps] + copy_all_cps[mo_ca[k]+1:] #replacing the old point with the new one
            all_cps = copy.deepcopy(copy_all_cps)
            fiber_moves += 1
            #print mo_ca[k], 'A fiber moved'#, cc_all_cps
            #sum_close_dist[mo_ca] = 0.0

            #second =



    all_cps = random_points(all_cps)
    #sum_close_dist, all_close_centers = closest_cents(all_cps)
    
    '''resin rich pocket'''
    resin_center = [h/2, w/2]
#    resin_rad = 2.0
#    res_dist = (resin_center, all_cps)
#    if min(res_dist) > resin_rad:
#            copy_all_cps = copy_all_cps[:mo_ca[k]] + [nps] + copy_all_cps[mo_ca[k]+1:] #replacing the old point with the new one
#            all_cps = copy.deepcopy(copy_all_cps)
    
    
    cp = len(all_cps)
    print cp
    add_step += 1
    if cp >= int(max_cp):
        break


#vf = cp*pi*r**2/(h*w)
#print vf









