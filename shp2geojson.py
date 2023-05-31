#!/usr/bin/env python
# coding: utf-8

__author__ = "Lasith Niroshan"
__copyright__ = "Copyright 2022, DeepMapper"
__credits__ = ["James D. Carswell"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Lasith Niroshan"
__email__ = "lasith.lan@gmail.com"
__status__ = "Research-And-Development"

"""
How to run:
                python .\shp2geojson.py --shapefile .\data\Export_Output.shp --db .\osm_database\ --feature buildings
"""

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import utm
import random
from shapely import geometry
import os, sys
import pandas, csv, json
from tqdm import tqdm

from math import radians, sqrt, sin, tan, cos, degrees


def arcmer(a, equad, lat1, lat2):
    b = a * sqrt(1 - equad)
    n = (a - b) / (a + b)
    a0 = 1.0 + ((n**2) / 4.0) + ((n**4) / 64.0)
    a2 = (3.0 / 2.0) * (n - ((n**3) / 8.0))
    a4 = (15.0 / 16.0) * ((n**2) - ((n**4) / 4.0))
    a6 = (35.0 / 48.0) * (n**3)
    s1 = (
        a
        / (1 + n)
        * (
            a0 * lat1
            - a2 * sin(2.0 * lat1)
            + a4 * sin(4.0 * lat1)
            - a6 * sin(6.0 * lat1)
        )
    )
    s2 = (
        a
        / (1 + n)
        * (
            a0 * lat2
            - a2 * sin(2.0 * lat2)
            + a4 * sin(4.0 * lat2)
            - a6 * sin(6.0 * lat2)
        )
    )
    return s2 - s1


def xy2geo(m, p, a, equad, lat0, lon0):
    lat0 = radians(lat0)
    lon0 = radians(lon0)
    sigma1 = p
    fil = lat0 + sigma1 / (a * (1 - equad))
    deltafi = 1
    while deltafi > 0.0000000001:
        sigma2 = arcmer(a, equad, lat0, fil)
        RO = a * (1 - equad) / ((1 - equad * (sin(fil) ** 2)) ** (3.0 / 2.0))
        deltafi = (sigma1 - sigma2) / RO
        fil = fil + deltafi
    N = a / sqrt(1 - equad * (sin(fil)) ** 2)
    RO = a * (1 - equad) / ((1 - equad * (sin(fil) ** 2)) ** (3.0 / 2.0))
    t = tan(fil)
    psi = N / RO
    lat = (
        fil
        - (t / RO) * ((m**2) / (2.0 * N))
        + (t / RO)
        * ((m**4) / (24.0 * (N**3)))
        * (-4.0 * (psi**2) - 9.0 * psi * (1.0 - t**2) + 12.0 * (t**2))
        - (t / RO)
        * (m**6 / (720.0 * (N**5)))
        * (
            8.0 * (psi**4) * (11.0 - 24.0 * (t**2))
            - 12.0 * (psi**3) * (21.0 - 71.0 * (t**2))
            + 15.0 * (psi**2) * (15.0 - 98.0 * (t**2) + 15.0 * (t**4))
            + 180.0 * psi * (5.0 * (t**2) - 3.0 * (t**4))
            - 360.0 * (t**4)
        )
        + (t / RO)
        * ((m**8) / (40320.0 * (N**7)))
        * (1385.0 + 3633.0 * (t**2) + 4095.0 * (t**4) + 1575.0 * (t**6))
    )
    lon = (
        (m / (N))
        - ((m**3) / (6.0 * (N**3))) * (psi + 2.0 * (t**2))
        + ((m**5) / (120.0 * (N**5)))
        * (
            -4.0 * (psi**3) * (1.0 - 6.0 * (t**2))
            + (psi**2) * (9.0 - 68.0 * (t**2))
            + 72.0 * psi * (t**2)
            + 24.0 * (t**4)
        )
        - ((m**7) / (5040.0 * (N**7)))
        * (61.0 + 662.0 * (t**2) + 1320.0 * (t**4) + 720.0 * (t**6))
    )
    lon = lon0 + lon / cos(fil)
    lat = degrees(lat)
    lon = degrees(lon)
    return lat, lon


def itm2geo(x, y):
    # GRS-80
    a = 6378137.0
    equad = 0.00669437999
    # Natural Origin
    lat0 = 53.5
    lon0 = -8.0
    k0 = 0.999820
    p = (y - 750000.0) / k0
    m = (x - 600000.0) / k0
    lat, lon = xy2geo(m, p, a, equad, lat0, lon0)
    return lat, lon


def load_map(path):
    """Load a geojson file.
    Args:
        path (str): path to the geojson file.
    Returns:
        object: geopandas object
    """
    geo_data = gpd.read_file(path)
    return geo_data


def buildGeoJson(arr, osi_id):
    geoJSON = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {"type": "Polygon", "coordinates": [arr]},
                "properties": {"provider": "OSI(c) 2021", "ID": osi_id},
            }
        ],
    }
    return geoJSON


def shp2geojsons(osm_database, feature, geoDF, ind):
    # coords = np.array(geoDF.iloc[ind]["geometry"].exterior).tolist()
    coords = np.array(geoDF.iloc[ind]["geometry"]).tolist()
    poly = [itm2geo(x[0], x[1]) for x in coords]
    poly2 = [(x[1], x[0]) for x in poly]
    geoObject = buildGeoJson(poly2, ind)
    # print(json.dumps(geoObject, indent=4))

    with open(f"{osm_database}/{feature}/{ind}.geojson", "w", encoding="utf-8") as f:
        json.dump(geoObject, f, ensure_ascii=False, indent=4)
    f.close()


def geojson2csv(in_file):
    out_file = in_file.replace("geojson", "csv")
    with open(in_file) as data:
        geojson = json.load(data)
        coords = geojson["features"][0]["geometry"]["coordinates"][0]
        with open(out_file, "wt") as csv:
            for node in coords:
                csv.write("%0.7f\t%0.7f\n" % (node[1], node[0]))
        csv.close()
    data.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Converts (ArcGIS) shape file into a collection of geojsons."
    )
    parser.add_argument("--shapefile", help="path for the shape file.")
    parser.add_argument("--db", help="path for the osm database.")
    parser.add_argument("--feature", help="feature name for the geojson.")
    args = parser.parse_args()
    shp_file = args.shapefile
    db_root = args.db
    feature = args.feature
    geoDF = load_map(shp_file)

    if not (os.path.exists(os.path.join(db_root, feature))):
        os.mkdir(os.path.join(db_root, feature))
    else:
        print(f"{db_root}/{feature} already exist!/nExiting...")
        exit(0)

    print("Converting shape file into geojsons....")

    for i in tqdm(range(0, geoDF.shape[0])):
        try:
            shp2geojsons(db_root, feature, geoDF, i)
        except Exception as expt:
            print("Error: ", expt)

    print("Building csv files for nodes.... ")

    db_dir = f"{db_root}/{feature}"
    geojsons = os.listdir(db_dir)
    for f in tqdm(geojsons):
        infile = f"{db_dir}/{f}"
        geojson2csv(infile)


if __name__ == "__main__":
    main()
