from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from caeModules import *
from driverUtils import executeOnCaeStartup

import random
import numpy as np
import csv
import copy
#import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt

session.journalOptions.setValues(replayGeometry=COORDINATE)

#set work directory
directory = "D:/1_res_paper/damage/2D/resin compare_3rd/0trial"



pl_width = 20.0
pl_length = 20.0
w = pl_width
h = pl_length
vf = 0.60
r = h/40.0
del_min = 0.01*r
seed_size = r/3
strain = 0.1  #percent
deformation = -1*strain/100.0 * w
#Carbon fiber transverse direction (MPa)
Ef, nuf = 28000.0, 0.23
Em, num = 3900.0, 0.34 #(epoxy) from fiedler

Sc = 114.5 #MPa
St = 47.0  #MPa

work_directory_n = 'D:/1_res_paper/damage/2D/resin compare_3rd/0trial/no res'
work_directory_w = 'D:/1_res_paper/damage/2D/resin compare_3rd/0trial/with res'



number_cps_jobs = []
rand_numb = 0
while rand_numb < 1:

    os.chdir(work_directory_n)
    ''' create part '''
    model_name = 'Model-1'
    part_name = 'Part-1'
    s = mdb.models[model_name].ConstrainedSketch(name='__profile__',
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(-pl_width/2, -pl_length/2), point2=(pl_width/2, pl_length/2))
    p = mdb.models[model_name].Part(name=part_name, dimensionality=TWO_D_PLANAR,
        type=DEFORMABLE_BODY)
    p = mdb.models[model_name].parts[part_name]
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()
    del mdb.models[model_name].sketches['__profile__']

    '''finding center points'''
    fiber_area = h*w*vf
    single_fiber = pi*r**2
    max_cp = round(fiber_area/single_fiber)
    execfile(directory+"/center_points.py")
    if cp != max_cp: # to make sure that Vf is exactly defined
        continue


    ''' making the partitions'''
    p = mdb.models[model_name].parts[part_name]
    f, e, d = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
        0.0, 0.0, 0.0))
    s = mdb.models[model_name].ConstrainedSketch(name='__profile__',
        sheetSize=150.0, gridSpacing=4.0, transform=t)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models[model_name].parts[part_name]
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

    for combined_points in all_cps:
        for single in combined_points:
            s.CircleByCenterPerimeter(center=single, point1=(single[0]+r, single[1]))
            
    
    points_del = []
    for combined_points in all_cps:
        for single in combined_points:
            s.CircleByCenterPerimeter(center=single, point1=(single[0]+r, single[1]))
            points_del.append(single)


    f = p.faces
    pickedFaces = f
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del mdb.models[model_name].sketches['__profile__']

    ''' material properties '''
    mdb.models[model_name].Material(name='Fiber')
    mdb.models[model_name].materials['Fiber'].Elastic(table=((Ef, nuf), ))
    mdb.models[model_name].HomogeneousSolidSection(name='Section-Fiber',
        material='Fiber', thickness=None)

    mdb.models[model_name].Material(name='Matrix')
    mdb.models[model_name].materials['Matrix'].Elastic(table=((Em, num), ))
    mdb.models[model_name].HomogeneousSolidSection(name='Section-Matrix',
        material='Matrix', thickness=None)
    mdb.models[model_name].materials['Matrix'].UserOutputVariables(n=1)

    # assigning sections
    p = mdb.models[model_name].parts[part_name]
    faces = p.faces
    region = regionToolset.Region(faces=faces)
    p.SectionAssignment(region=region, sectionName='Section-Matrix', offset=0.0,
        offsetType=MIDDLE_SURFACE, offsetField='',
        thicknessAssignment=FROM_SECTION)

    #fiber section assignment
    p = mdb.models[model_name].parts[part_name]
    f = p.faces
    eps = r*1e-5
    for center_points in all_cps:
        for cpo in center_points:
            #if cp[0] > -w/2 and cp[0] < w/2 and cp[1] > -h/2 and cp[1] < h/2:
            xf, yf = cpo[0], cpo[1]
            if cpo[0] < -w/2:
                xf = -w/2+eps
            if cpo[0] > w/2:
                xf = w/2-eps
            if cpo[1] < -h/2:
                yf = -h/2+eps
            if cpo[1] > h/2:
                yf = h/2-eps

            faces = f.findAt(((xf,yf,0.0),))
            region = regionToolset.Region(faces=faces)
            p.SectionAssignment(region=region, sectionName='Section-Fiber', offset=0.0,
                offsetType=MIDDLE_SURFACE, offsetField='',
                thicknessAssignment=FROM_SECTION)

    ''' Mesh '''
    p = mdb.models[model_name].parts[part_name]
    f = p.faces
    pickedRegions = f
    
    elemType1 = mesh.ElemType(elemCode=CPE4R, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, hourglassControl=DEFAULT, 
        distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT)
#p = mdb.models['Model-1'].parts['Part-1']
#f = p.faces
#pickedRegions =(faces, )
#p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
#p = mdb.models['Model-1'].parts['Part-1']
#p.seedPart(size=0.24, deviationFactor=0.1, minSizeFactor=0.1)
#p = mdb.models['Model-1'].parts['Part-1']
#p.generateMesh()   
    pickedRegion2 = (f, )
    p.setElementType(regions=(f, ), elemTypes=(elemType1, elemType2))
    p.setMeshControls(regions=pickedRegions, elemShape=TRI)
    p = mdb.models[model_name].parts[part_name]
    p.seedPart(size=seed_size, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models[model_name].parts[part_name]
    p.generateMesh()

    #=====================
    ''' creating step '''
    mdb.models[model_name].StaticStep(name='Step-1', previous='Initial')
    mdb.models[model_name].fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'E', 'UVARM', 'U'))

    ''' assembly '''
    a = mdb.models[model_name].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models[model_name].parts[part_name]
    insta_name = part_name + '-1'
    a.Instance(name=insta_name, part=p, dependent=ON)

    ''' Periodic Boundary Conditions'''
    p = mdb.models[model_name].parts[part_name]
    #execfile(directory+"\\periodic.py")
    #%%%% this part is for PBC
    #dummy nodes for applying load and BC
    a = mdb.models[model_name].rootAssembly
    refp_load = (w/2+1.0, 0.0, 0.0)
    refp_bc = (0.0, h/2+1.0, 0.0)
    
    a.ReferencePoint(point=refp_load)
    a.ReferencePoint(point=refp_bc)
    
    
    #a = mdb.models['Model-1'].rootAssembly
    rr = mdb.models[model_name].rootAssembly.referencePoints
    refPoints=(rr[4], )
    a.Set(referencePoints=refPoints, name='Set-Ref-Load')
    refPoints=(rr[5], )
    a.Set(referencePoints=refPoints, name='Set-Ref-BC')
    
    #########
    left_nodes  = p.nodes.getByBoundingBox(-w/2-eps,-h/2-eps,0.0,  -w/2+eps, h/2+eps,0.0)
    right_nodes = p.nodes.getByBoundingBox(w/2-eps,-h/2-eps,0.0,  w/2+eps, h/2+eps,0.0)
    top_nodes   = p.nodes.getByBoundingBox(-w/2-eps, h/2-eps,0.0,  w/2+eps, h/2+eps,0.0)
    bot_nodes   = p.nodes.getByBoundingBox(-w/2-eps, -h/2-eps,0.0,  w/2+eps, -h/2+eps,0.0)
    
    
    if len(left_nodes) != len(right_nodes) or len(top_nodes) != len(bot_nodes): # for those cases that number of nodes aren't exactly equal
        Mdb()
        continue
          
    #left node coordinates
    left_node_co_y = []
    for i in range(len(left_nodes)):
        left_node_co_y.append(left_nodes[i].coordinates[1])
    
    left_nodes_org = []
    ll = len(left_node_co_y)
    for i in range(ll):
        min_index = left_node_co_y.index(min(left_node_co_y))
        left_nodes_org.append(left_nodes[min_index].label)
        left_node_co_y[min_index] = []
    
    
    #right node coordinates
    right_node_co_y = []
    for i in range(len(right_nodes)):
        right_node_co_y.append(right_nodes[i].coordinates[1])
    
    right_nodes_org = []
    ll = len(right_node_co_y)
    for i in range(ll):
        min_index = right_node_co_y.index(min(right_node_co_y))
        right_nodes_org.append(right_nodes[min_index].label)
        right_node_co_y[min_index] = []
    
    
    #top node coordinates
    top_node_co_x = []
    for i in range(len(top_nodes)):
        top_node_co_x.append(top_nodes[i].coordinates[0])
    
    top_nodes_org = []
    ll = len(top_node_co_x)
    for i in range(ll):
        min_index = top_node_co_x.index(min(top_node_co_x))
        top_nodes_org.append(top_nodes[min_index].label)
        top_node_co_x[min_index] = []
    
    
    #bot node coordinates
    bot_node_co_x = []
    for i in range(len(bot_nodes)):
        bot_node_co_x.append(bot_nodes[i].coordinates[0])
    
    bot_nodes_org = []
    ll = len(bot_node_co_x)
    for i in range(ll):
        min_index = bot_node_co_x.index(min(bot_node_co_x))
        bot_nodes_org.append(bot_nodes[min_index].label)
        bot_node_co_x[min_index] = []
    
    
    #equation boundary conditions
    #create set
    a = mdb.models[model_name].rootAssembly
    a.regenerate()
    n = a.instances[insta_name].nodes
    for i in range(len(left_nodes_org)):
        a.Set(nodes=n[left_nodes_org[i]-1:left_nodes_org[i]], name='set-left-'+str(i))
        a.Set(nodes=n[right_nodes_org[i]-1:right_nodes_org[i]], name='set-right-'+str(i))
    
    for i in range(len(bot_nodes_org)):
        a.Set(nodes=n[bot_nodes_org[i]-1:bot_nodes_org[i]], name='set-bot-'+str(i))
        a.Set(nodes=n[top_nodes_org[i]-1:top_nodes_org[i]], name='set-top-'+str(i))
    
    
    for dim in [1,2]:
        for i in range(len(left_nodes_org)):
            mdb.models[model_name].Equation(name='leftright-dim'+str(dim)+'-node-'+str(i), terms=((-1.0, 'set-left-'+str(i), dim),
                (1.0, 'set-right-'+str(i), dim), (1.0, 'Set-Ref-Load', dim)))
    
    for dim in [1,2]:
        for i in range(len(top_nodes_org)-2):
            mdb.models[model_name].Equation(name='topbot-dim'+str(dim)+'-node-'+str(i+1), terms=((-1.0, 'set-bot-'+str(i+1), dim),
                (1.0, 'set-top-'+str(i+1), dim), (1.0, 'Set-Ref-BC', dim)))
    
    
    #%%%
    
    ''' Load and BC '''
    a = mdb.models[model_name].rootAssembly
    region = a.sets['Set-Ref-BC']
    mdb.models[model_name].DisplacementBC(name='BC-dis', createStepName='Step-1',
        region=region, u1=0.0, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF,
        distributionType=UNIFORM, fieldName='', localCsys=None)

    region = a.sets['Set-Ref-Load']
    mdb.models[model_name].DisplacementBC(name='BC-Load', createStepName='Step-1',
        region=region, u1=deformation, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF,
        distributionType=UNIFORM, fieldName='', localCsys=None)


    for strain in np.arange(0.1,0.4,0.2): #applying 0.1 and 0.3 strain
        deformation = -1*strain/100.0 * w
        mdb.models['Model-1'].boundaryConditions['BC-Load'].setValues(u1=deformation)       
                
        '''submitting job '''
        job_name='Random-' + str(rand_numb) + '-strain-' + str(int(strain*1000))
        mdb.Job(name=job_name, model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
        userSubroutine=directory+'\UVARM.for',
        scratch='', parallelizationMethodExplicit=DOMAIN,
        multiprocessingMode=DEFAULT, numDomains=1, numCpus=1)
        mdb.jobs[job_name].submit(consistencyChecking=OFF)
        mdb.jobs[job_name].waitForCompletion()
    
        ''' Post production '''
        #number_cps_jobs.append([rand_numb, cp])
        execfile(directory+"/postproduction_max_par.py")


    fiber_loc_name = '/fiber_loc_' +str(rand_numb)+ '.csv'
    file_name = newpath +  fiber_loc_name
    with open( file_name, 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        #wr.writerow(('mat area','mat S11', 'mat E11'))#, 'mat S12'))
        for row in points_del:    
            wr.writerow(row)            
        
    
    
    Mdb()
# doing samples with resin pocket  
#''' another two analysis'''
#############################################################################
####
####
#############################################################################    
#############################################################################
####
####
#############################################################################
#############################################################################
####
####
#############################################################################
    os.chdir(work_directory_w)
    ''' create part '''
    model_name = 'Model-1'
    part_name = 'Part-2'
    s = mdb.models[model_name].ConstrainedSketch(name='__profile__',
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(-pl_width/2, -pl_length/2), point2=(pl_width/2, pl_length/2))
    p = mdb.models[model_name].Part(name=part_name, dimensionality=TWO_D_PLANAR,
        type=DEFORMABLE_BODY)
    p = mdb.models[model_name].parts[part_name]
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()
    del mdb.models[model_name].sketches['__profile__']

    '''finding center points'''
    fiber_area = h*w*vf
    single_fiber = pi*r**2
    max_cp = round(fiber_area/single_fiber)
    execfile(directory+"/center_points_elim.py")
#    if cp != max_cp: # to make sure that Vf is exactly defined
#        continue


    ''' making the partitions'''
    p = mdb.models[model_name].parts[part_name]
    f, e, d = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
        0.0, 0.0, 0.0))
    s = mdb.models[model_name].ConstrainedSketch(name='__profile__',
        sheetSize=150.0, gridSpacing=4.0, transform=t)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models[model_name].parts[part_name]
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

    for combined_points in all_cps:
        for single in combined_points:
            s.CircleByCenterPerimeter(center=single, point1=(single[0]+r, single[1]))
            
    
    points_del = []
    for combined_points in all_cps:
        for single in combined_points:
            s.CircleByCenterPerimeter(center=single, point1=(single[0]+r, single[1]))
            points_del.append(single)


    f = p.faces
    pickedFaces = f
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del mdb.models[model_name].sketches['__profile__']

    ''' material properties '''
    mdb.models[model_name].Material(name='Fiber')
    mdb.models[model_name].materials['Fiber'].Elastic(table=((Ef, nuf), ))
    mdb.models[model_name].HomogeneousSolidSection(name='Section-Fiber',
        material='Fiber', thickness=None)

    mdb.models[model_name].Material(name='Matrix')
    mdb.models[model_name].materials['Matrix'].Elastic(table=((Em, num), ))
    mdb.models[model_name].HomogeneousSolidSection(name='Section-Matrix',
        material='Matrix', thickness=None)
    mdb.models[model_name].materials['Matrix'].UserOutputVariables(n=1)

    # assigning sections
    p = mdb.models[model_name].parts[part_name]
    faces = p.faces
    region = regionToolset.Region(faces=faces)
    p.SectionAssignment(region=region, sectionName='Section-Matrix', offset=0.0,
        offsetType=MIDDLE_SURFACE, offsetField='',
        thicknessAssignment=FROM_SECTION)

    #fiber section assignment
    p = mdb.models[model_name].parts[part_name]
    f = p.faces
    eps = r*1e-5
    for center_points in all_cps:
        for cpo in center_points:
            #if cp[0] > -w/2 and cp[0] < w/2 and cp[1] > -h/2 and cp[1] < h/2:
            xf, yf = cpo[0], cpo[1]
            if cpo[0] < -w/2:
                xf = -w/2+eps
            if cpo[0] > w/2:
                xf = w/2-eps
            if cpo[1] < -h/2:
                yf = -h/2+eps
            if cpo[1] > h/2:
                yf = h/2-eps

            faces = f.findAt(((xf,yf,0.0),))
            region = regionToolset.Region(faces=faces)
            p.SectionAssignment(region=region, sectionName='Section-Fiber', offset=0.0,
                offsetType=MIDDLE_SURFACE, offsetField='',
                thicknessAssignment=FROM_SECTION)

    ''' Mesh '''
    p = mdb.models[model_name].parts[part_name]
    f = p.faces
    pickedRegions = f
    
    elemType1 = mesh.ElemType(elemCode=CPE4R, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, hourglassControl=DEFAULT, 
        distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT)
#p = mdb.models['Model-1'].parts['Part-1']
#f = p.faces
#pickedRegions =(faces, )
#p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
#p = mdb.models['Model-1'].parts['Part-1']
#p.seedPart(size=0.24, deviationFactor=0.1, minSizeFactor=0.1)
#p = mdb.models['Model-1'].parts['Part-1']
#p.generateMesh()   
    pickedRegion2 = (f, )
    p.setElementType(regions=(f, ), elemTypes=(elemType1, elemType2))
    p.setMeshControls(regions=pickedRegions, elemShape=TRI)
    p = mdb.models[model_name].parts[part_name]
    p.seedPart(size=seed_size, deviationFactor=0.1, minSizeFactor=0.1)
    p = mdb.models[model_name].parts[part_name]
    p.generateMesh()

    #=====================
    ''' creating step '''
    mdb.models[model_name].StaticStep(name='Step-1', previous='Initial')
    mdb.models[model_name].fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'E', 'UVARM', 'U'))

    ''' assembly '''
    a = mdb.models[model_name].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models[model_name].parts[part_name]
    insta_name = part_name + '-1'
    a.Instance(name=insta_name, part=p, dependent=ON)

    ''' Periodic Boundary Conditions'''
    p = mdb.models[model_name].parts[part_name]
    #execfile(directory+"\\periodic.py")
    #%%%% this part is for PBC
    #dummy nodes for applying load and BC
    a = mdb.models[model_name].rootAssembly
    refp_load = (w/2+1.0, 0.0, 0.0)
    refp_bc = (0.0, h/2+1.0, 0.0)
    
    a.ReferencePoint(point=refp_load)
    a.ReferencePoint(point=refp_bc)
    
    
    #a = mdb.models['Model-1'].rootAssembly
    rr = mdb.models[model_name].rootAssembly.referencePoints
    refPoints=(rr[4], )
    a.Set(referencePoints=refPoints, name='Set-Ref-Load')
    refPoints=(rr[5], )
    a.Set(referencePoints=refPoints, name='Set-Ref-BC')
    
    #########
    left_nodes  = p.nodes.getByBoundingBox(-w/2-eps,-h/2-eps,0.0,  -w/2+eps, h/2+eps,0.0)
    right_nodes = p.nodes.getByBoundingBox(w/2-eps,-h/2-eps,0.0,  w/2+eps, h/2+eps,0.0)
    top_nodes   = p.nodes.getByBoundingBox(-w/2-eps, h/2-eps,0.0,  w/2+eps, h/2+eps,0.0)
    bot_nodes   = p.nodes.getByBoundingBox(-w/2-eps, -h/2-eps,0.0,  w/2+eps, -h/2+eps,0.0)
    
    
    if len(left_nodes) != len(right_nodes) or len(top_nodes) != len(bot_nodes): # for those cases that number of nodes aren't exactly equal
        Mdb()
        continue
          
    #left node coordinates
    left_node_co_y = []
    for i in range(len(left_nodes)):
        left_node_co_y.append(left_nodes[i].coordinates[1])
    
    left_nodes_org = []
    ll = len(left_node_co_y)
    for i in range(ll):
        min_index = left_node_co_y.index(min(left_node_co_y))
        left_nodes_org.append(left_nodes[min_index].label)
        left_node_co_y[min_index] = []
    
    
    #right node coordinates
    right_node_co_y = []
    for i in range(len(right_nodes)):
        right_node_co_y.append(right_nodes[i].coordinates[1])
    
    right_nodes_org = []
    ll = len(right_node_co_y)
    for i in range(ll):
        min_index = right_node_co_y.index(min(right_node_co_y))
        right_nodes_org.append(right_nodes[min_index].label)
        right_node_co_y[min_index] = []
    
    
    #top node coordinates
    top_node_co_x = []
    for i in range(len(top_nodes)):
        top_node_co_x.append(top_nodes[i].coordinates[0])
    
    top_nodes_org = []
    ll = len(top_node_co_x)
    for i in range(ll):
        min_index = top_node_co_x.index(min(top_node_co_x))
        top_nodes_org.append(top_nodes[min_index].label)
        top_node_co_x[min_index] = []
    
    
    #bot node coordinates
    bot_node_co_x = []
    for i in range(len(bot_nodes)):
        bot_node_co_x.append(bot_nodes[i].coordinates[0])
    
    bot_nodes_org = []
    ll = len(bot_node_co_x)
    for i in range(ll):
        min_index = bot_node_co_x.index(min(bot_node_co_x))
        bot_nodes_org.append(bot_nodes[min_index].label)
        bot_node_co_x[min_index] = []
    
    
    #equation boundary conditions
    #create set
    a = mdb.models[model_name].rootAssembly
    a.regenerate()
    n = a.instances[insta_name].nodes
    for i in range(len(left_nodes_org)):
        a.Set(nodes=n[left_nodes_org[i]-1:left_nodes_org[i]], name='set-left-'+str(i))
        a.Set(nodes=n[right_nodes_org[i]-1:right_nodes_org[i]], name='set-right-'+str(i))
    
    for i in range(len(bot_nodes_org)):
        a.Set(nodes=n[bot_nodes_org[i]-1:bot_nodes_org[i]], name='set-bot-'+str(i))
        a.Set(nodes=n[top_nodes_org[i]-1:top_nodes_org[i]], name='set-top-'+str(i))
    
    
    for dim in [1,2]:
        for i in range(len(left_nodes_org)):
            mdb.models[model_name].Equation(name='leftright-dim'+str(dim)+'-node-'+str(i), terms=((-1.0, 'set-left-'+str(i), dim),
                (1.0, 'set-right-'+str(i), dim), (1.0, 'Set-Ref-Load', dim)))
    
    for dim in [1,2]:
        for i in range(len(top_nodes_org)-2):
            mdb.models[model_name].Equation(name='topbot-dim'+str(dim)+'-node-'+str(i+1), terms=((-1.0, 'set-bot-'+str(i+1), dim),
                (1.0, 'set-top-'+str(i+1), dim), (1.0, 'Set-Ref-BC', dim)))
    
    
    #%%%
    
    ''' Load and BC '''
    a = mdb.models[model_name].rootAssembly
    region = a.sets['Set-Ref-BC']
    mdb.models[model_name].DisplacementBC(name='BC-dis', createStepName='Step-1',
        region=region, u1=0.0, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF,
        distributionType=UNIFORM, fieldName='', localCsys=None)

    region = a.sets['Set-Ref-Load']
    mdb.models[model_name].DisplacementBC(name='BC-Load', createStepName='Step-1',
        region=region, u1=deformation, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF,
        distributionType=UNIFORM, fieldName='', localCsys=None)


    for strain in np.arange(0.1,0.4,0.2): #applying 0.1 and 0.3 strain
        deformation = -1*strain/100.0 * w
        mdb.models['Model-1'].boundaryConditions['BC-Load'].setValues(u1=deformation)       
                
        '''submitting job '''
        job_name='Random-res' + str(rand_numb) + '-strain-' + str(int(strain*1000))
        mdb.Job(name=job_name, model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
        userSubroutine=directory+'\UVARM.for',
        scratch='', parallelizationMethodExplicit=DOMAIN,
        multiprocessingMode=DEFAULT, numDomains=1, numCpus=1)
        mdb.jobs[job_name].submit(consistencyChecking=OFF)
        mdb.jobs[job_name].waitForCompletion()
    
        ''' Post production '''
        number_cps_jobs.append([rand_numb, len(all_cps)])
        execfile(directory+"/postproduction_max_par_elim.py")


    fiber_loc_name = '/fiber_loc_' +str(rand_numb)+ '.csv'
    file_name = newpath +  fiber_loc_name
    with open( file_name, 'wb') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        #wr.writerow(('mat area','mat S11', 'mat E11'))#, 'mat S12'))
        for row in points_del:    
            wr.writerow(row)            
 

    
    rand_numb += 1
    Mdb()#delete the model data base
    


file_name = newpath +  'all_number_of_fibers.csv'
with open( file_name, 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for row in number_cps_jobs:    
        wr.writerow(row) 











