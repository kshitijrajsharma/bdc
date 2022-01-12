#!/usr/bin/env python
# coding: utf-8

# In[1]:


import geopandas as gpd
import pandas as pd
import pandas
import math
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon
from geographiclib.geodesic import Geodesic
from shapely.ops import split
from shapely import ops
from pyproj import Proj, Transformer


# In[2]:


buildings= gpd.read_file('../tests/data/buildings.geojson')
road= gpd.read_file('../tests/data/road.geojson')


# In[3]:


# for d in buildings.columns:
#     print(d)
road['end'] = None
road['start'] = None
road['start_point_x'] = None
road['start_point_y'] = None
road['end_point_x'] = None
road['end_point_y'] = None
road['length'] = None
road['bearing'] = None
road['intersected'] = None
road['mother_road']='None'

buildings['road_direction'] = 'None'
buildings['distance2_intersecting_road'] = 0
buildings['road_length_d_code'] = ''

# print(buildings)
# print(road)


# In[4]:


class bdcPoint:
    	
	def __init__(self):
		
		self.x = 0
		self.y = 0

# Constant integers for directions
RIGHT = 1
LEFT = -1
ZERO = 0

def directionOfPoint(A, B, P):
    """[Provides direction of the point in respect to line]

    Args:
        A ([type]): [Starting Point of Line]
        B ([type]): [Ending Point of Line]
        P ([type]): [Point of Interest]

    Returns:
        [Direction(int)]: [Left,Right and Zero : -1,1,0 ( Left direction, Right Direction , On Line)]
    """
    global RIGHT, LEFT, ZERO
	
	# Subtracting co-ordinates of
	# point A from B and P, to
	# make A as origin
 
    B.x -= A.x
    B.y -= A.y
    P.x -= A.x
    P.y -= A.y

    # Determining cross Product
    """"The Cross-Product has an interesting Property which will be used to determine direction of a point from a line segment. That is, the cross-product of two points is positive if and only if the angle of those point at origin (0, 0) is in counter-clockwise. And conversely the cross-product is negative if and only if the angle of those point at origin is in clockwise direction."""
    cross_product = B.x * P.y - B.y * P.x
    # Return RIGHT if cross product is positive
    if (cross_product > 0):
        return RIGHT
        
    # Return LEFT if cross product is negative
    if (cross_product < 0):
        return LEFT

    # Return ZERO if cross product is zero
    return ZERO

def addDirection_BuildingCode(direction,building_length,road_length,index,road_index):
    #right-even, left-odd number 
    # print(length)
    building_length=int(building_length)
    if road_length < 60:
        associated_road=road.iloc[road_index]

        if associated_road.mother_road == 'start':
            mother_road_code=associated_road.start
        else :
            mother_road_code=associated_road.end

        mother_road= road.loc[road['NEW'] == mother_road_code]

        mr_geom= mother_road.geometry
        d = {'col1': ['name1'], 'geometry': [Point(associated_road.start_point_x, associated_road.start_point_y)]}
        motherroad_start_point = {'col1': ['name1'], 'geometry': [Point(mother_road.start_point_x, mother_road.start_point_y)]}
        
        gdf = gpd.GeoDataFrame(d, crs=4326)
        mother_df = gpd.GeoDataFrame(motherroad_start_point, crs=4326)
        
        sp_point=gdf.iloc[0]
        mother_st_p=mother_df.iloc[0]
        for line in mr_geom:
            # print(line)
            for l in line:
                mother_line=l
                break

        result=split(mother_line,sp_point.geometry)
        for r  in result:
            splitted_line=r
            splitted_line.srid = 4326
            if splitted_line.contains(mother_st_p.geometry):
                break
            else:
                continue
        my_transformer = Transformer.from_crs('EPSG:4326', 'EPSG:32644', always_xy=True)
        geom_transformed = ops.transform(my_transformer.transform, splitted_line)
        
        buildings.at[index,'road_length_d_code']=str(int(geom_transformed.length))+"/"
        if (direction == 1):
            # print("Right Direction")
            buildings.at[index,'road_direction']="right"
            buildings.at[index,'road_length_d_code']+=str(int(round_up_to_even(building_length)))

        elif (direction == -1):
            # print("Left Direction")
            buildings.at[index,'road_direction']="left"
            buildings.at[index,'road_length_d_code']+=str(int(round_up_to_odd(building_length)))
        else:
            print("On  Line")
            buildings.at[index,'road_direction']="on_line"
    else :
        if (direction == 1):
            # print("Right Direction")
            buildings.at[index,'road_direction']="right"
            buildings.at[index,'road_length_d_code']=str(int(round_up_to_even(building_length)))

        elif (direction == -1):
            # print("Left Direction")
            buildings.at[index,'road_direction']="left"
            buildings.at[index,'road_length_d_code']=str(int(round_up_to_odd(building_length)))
        else:
            print("On  Line")
            buildings.at[index,'road_direction']="on_line"
    
    
#     print(buildings.at[index,'road_length_d_code'])
            
def getStartEndPoint(line,type):
    if len(line)==1 or type=="road" :
       
        df_line = line.geometry.boundary
        # print(len(df_line))
        # print(line_bound)
        # df_line=line_bound.explode(index_parts=True)
        start_point,end_point=df_line[0],df_line[1]

        return start_point,end_point
    else:
        raise ValueError ("Line contains more than one code")

def getLineLength(line):
    projected_line = line.to_crs(epsg=32644)
    length= projected_line.geometry.length
    float_length=length.tolist()
    return float_length

def round_up_to_even(f):
    return math.ceil(f / 2.) * 2

def round_up_to_odd(f):
    return math.ceil(f) // 2 * 2 + 1

def addStartEndPoint(start_point,end_point,index):
    
    road.at[index,'start_point_x']=float(start_point.x)
    road.at[index,'start_point_y']=float(start_point.y)
    
    road.at[index,'end_point_x']=float(end_point.x)
    road.at[index,'end_point_y']=float(end_point.y)

def ReviseStartEndDateWithBearing(bearing,index):
    # default direction is west to east , south to north 
    # print("ma revise garna aako dubai vettera")
    if bearing is None :
        bearing=0
    bearing=int(bearing)
    def assign_start():
        print("start")
        
        road.at[index,'mother_road']='start'

    def assign_end():
    # east_west/west_east direction , assign to east
        print("end")
        road.at[index,'mother_road']='end'

    
    if ((bearing >= 0) and (bearing<=45)) :
        #    North_South/south_north Direction, assign to south
        
        assign_start()
        
    if ((bearing > 45) and (bearing<135)) :
        
        # east_west/west_east direction , assign to east
        assign_end()
    
    if ((bearing >= 135) and (bearing<=225)) :
        #    North_South/south_north Direction, assign to south
        
        assign_start()
    if ((bearing > 225) and (bearing<315)) :
            
        # east_west/west_east direction , assign to east
        assign_end()
    if ((bearing >= 315) and (bearing<=360)) :
    #    North_South/south_north Direction, assign to south
        
        assign_start()      

def FindIntersectedFeatures(start_point,end_point,line_row,index_road):
    for index, row in road.iterrows():
        line=row.geometry        
        if str(row.NEW) != str(line_row.NEW):
            if line.distance(start_point) < 1e-8 :
                # print("start found "+ str(row.NEW))
                # print("mero start vetiyo guys ")
                road.at[index_road,'start']=str(row.NEW)
    
            if line.distance(end_point) < 1e-8:
                # print("end found "+ str(row.NEW))
                road.at[index_road,'end']=str(row.NEW)
                # print("mero end vetiyo guys ")


def scan_I_feature(line_row,index_road):
    
        if ((line_row.start!= None) and (line_row.end== None)) :
            road.at[index_road,'mother_road']='start'
        if ((line_row.start== None) and (line_row.end != None)) :
            road.at[index_road,'mother_road']='end'
    
        if ((line_row.start!= None) and (line_row.end != None)) :
        #both start and end exist which means line is in between roads 
            road.at[index_road,'mother_road']='None'
            ReviseStartEndDateWithBearing(line_row.bearing,index_road)

                    
def CalculateBearing(start_point,end_point,index):
    brng = Geodesic.WGS84.Inverse(start_point.y, start_point.x, end_point.y, end_point.x)['azi1']
    road.at[index,'bearing']=brng

def is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False    
    
def AssignIntersectedPoint(start_point_assumed,end_point_assumed,row,index):
#     print(row)
    if  row.mother_road :
        mother_road_direction = row.mother_road
        mother_road_code=None
        if mother_road_direction == 'start':
            mother_road_code=row.start
        elif mother_road_direction == 'end' :
            mother_road_code=row.end
        mother_road= road.loc[road['NEW'] == mother_road_code].geometry
        current_road = row.geometry
        mother_road_geom=mother_road
        for l in mother_road:
            for single in l:
                mother_road_geom=single
            break
#         print(mother_road_geom)
        for l in current_road:
            road_geom=l
            break
#         print(road_geom)
        start_point=mother_road_geom.intersection(road_geom)
        print(start_point)
        st=str(type(start_point))
#         print(st)
#         geoseries=f"""<class 'geopandas.geoseries.GeoSeries'>"""
#         print(geoseries)
    
        if isinstance(start_point,pandas.core.series.Series) :
            print(type(start_point))
            print('not counted')
            
            road.at[index,'intersected']='No'
        else:
            print(start_point)
            try:
                x=start_point.x
                y=start_point.y
                road.at[index,'start_point_x']=float(x)
                road.at[index,'start_point_y']=float(y)              
            except:
                try:
                    for p in start_point:
                        x=p.x
                        y=p.y
                        break
                    road.at[index,'start_point_x']=x
                    road.at[index,'start_point_y']=y
                except:
                    try:
                        df_line = start_point.boundary
                        start=df_line[0]
                        road.at[index,'start_point_x']=start.x
                        road.at[index,'start_point_y']=start.y
                    except:
                        print('i didnt assigned')
                        pass
#                     road.at[index,'start_point_x']=0
#                     road.at[index,'start_point_y']=0
 

        
def get_building_projected_road_length(road_linked_wrt_point,row):
    road_code = row.NEW
#     print(road_linked_wrt_point)
    linked_road_geometry = road_linked_wrt_point.geometry
    building_geometry = row.geometry
    degree_length_for_building=linked_road_geometry.project(building_geometry)
    # 5m=0.00005 in degree 
    building_point_projected_on_road=linked_road_geometry.interpolate(degree_length_for_building) #end_always 
    for r in building_point_projected_on_road:
        building_point_projected_on_road_single_geom=r
        break
    try:
        projected_building_on_road = {'col1': ['name1'], 'geometry': [Point(building_point_projected_on_road_single_geom.x, building_point_projected_on_road_single_geom.y)]}
    except:
        print(projected_building_on_road)
        raise("found error")

    road_link_wrt_point_x=road_linked_wrt_point.start_point_x
    road_link_wrt_point_y=road_linked_wrt_point.start_point_y
    try:
        road_start_point = {'col1': ['name1'], 'geometry': [Point(road_link_wrt_point_x,road_link_wrt_point_y)]}

        projected_building_on_road_gdf = gpd.GeoDataFrame(projected_building_on_road, crs=4326)


        road_start_df = gpd.GeoDataFrame(road_start_point, crs=4326)


        splitting_point=projected_building_on_road_gdf.iloc[0]
        start_point=road_start_df.iloc[0]
#         print(start_point.geometry)
        for line in linked_road_geometry:
            # print(line)
            for l in line:
                mother_line=l
                break
        result=split(mother_line,splitting_point.geometry.buffer(0.0000001))

        for r  in result:
                splitted_line=r
                splitted_line.srid = 4326
                if splitted_line.buffer(0.00001).contains(start_point.geometry.buffer(0.0000001)) == True:

                    break
                else:
                    continue


        
        my_transformer = Transformer.from_crs('EPSG:4326', 'EPSG:32644', always_xy=True)
        geom_transformed = ops.transform(my_transformer.transform, splitted_line)

        clipped_length=geom_transformed.length
        return clipped_length ,splitted_line ,building_point_projected_on_road_single_geom,start_point.geometry
    except:
        pass



    
                   


# In[5]:


def CalculateTouchLine(index,row):
    
    # print(row)
    start_point,end_point=getStartEndPoint(row,"road")
    addStartEndPoint(start_point,end_point,index)
    CalculateBearing(start_point,end_point,index)
    FindIntersectedFeatures(start_point,end_point,row,index)
    scan_I_feature(row,index)
    AssignIntersectedPoint(start_point,end_point,row,index)
    # break


# In[12]:


for index,row in road.iterrows():
    CalculateTouchLine(index,row)
    print(index)
print("Sucessfully added : startpoint, endpoint , Start touching Line , End Touching Line , bearing ")   


# In[75]:


def assign_building_code(index,row):  

    road_linked_wrt_point= road.loc[road['NEW'] == row.NEW]
#     print(road_linked_wrt_point)
    for l in road_linked_wrt_point.geometry:
        line_geom=l
        break

    my_transformer = Transformer.from_crs('EPSG:4326', 'EPSG:32644', always_xy=True)
    geom_transformed = ops.transform(my_transformer.transform, line_geom)

    road_length=geom_transformed.length

    clipped_length,splitted_mother_road_side,building_point_projected_on_road_single_geom,start_point=get_building_projected_road_length(road_linked_wrt_point,row)

    splitted_with_buffer_boi=split(splitted_mother_road_side,building_point_projected_on_road_single_geom.buffer(0.0001)) 

    for r  in splitted_with_buffer_boi:
        roi_line=r
        roi_line.srid = 4326
        if roi_line.buffer(0.00001).contains(start_point.buffer(0.0000001)) == False:
            break
        else:
            continue

    df_line = roi_line.boundary
    start_roi,end_roi=df_line[0],df_line[1]
#     print(start_roi)
#     print(end_roi)
#     print(building_point_projected_on_road_single_geom)
    if (round(start_roi.x,6) == round(building_point_projected_on_road_single_geom.x,4)) and (round(start_roi.y,4) == round(building_point_projected_on_road_single_geom.y,4)) :
        start_roi=end_roi
        end_roi=building_point_projected_on_road_single_geom
        print("start_end changed")

    road_index=road_linked_wrt_point.index[0]
#         road.at[road_index,'length']=length[0]       
    start = bdcPoint()
    end = bdcPoint()
    point = bdcPoint()
    poi = row.geometry

    start.x,start.y,end.x,end.y,point.x,point.y=start_roi.x,start_roi.y,end_roi.x,end_roi.y,poi.x,poi.y

    direction=directionOfPoint(start,end,point)
    # print(direction)
    # break
#     print(clipped_length)
#     print(direction)
    buildings.at[index,'distance2_intersecting_road']=clipped_length
    road.at[road_index,'length']=road_length
    
    
    addDirection_BuildingCode(direction,clipped_length,road_length,index,road_index)
    
        
#     print("added building code")


# In[5]:


for index,row in buildings.iterrows():
    try:
        assign_building_code(index,row)
        print(buildings.iloc[index].distance2_intersecting_road)
        print(index)
    except:
        continue


# In[ ]:


for index, row in buildings.iterrows():
    print(row)
    try:
        assign_building_code(index,row)
        print(index)
        print(buildings.iloc[index].road_length_d_code)
    except:
        continue
print("-------Sucess ------ added: road_length_d_code , length of road , direction of building ")
    


# In[ ]:


buildings.to_file("output_building_1.geojson", driver="GeoJSON")


# In[ ]:


print(buildings)


# In[ ]:


road.to_file("output_road_1.geojson", driver="GeoJSON")


# In[ ]:





# In[ ]:




