import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon

def main():
    input_geojson= gpd.read_file()
    a = Point(0, 0)
    b = Point(1, 0)
    distance= a.distance(b)
    print(distance)

if __name__ == "__main__":
    main()
