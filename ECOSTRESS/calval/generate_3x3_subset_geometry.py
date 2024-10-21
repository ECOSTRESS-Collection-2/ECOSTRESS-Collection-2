import rasters as rt

def generate_3x3_subset_geometry(geometry: rt.RasterGeometry, point: rt.Point):
    tower_row, tower_col = geometry.index_point(point.to_crs(geometry.crs))
    rows, cols = geometry.shape
    subset_3x3 = geometry[max(tower_row - 1, 0):min(tower_row + 2, rows - 1), max(tower_col - 1, 0):min(tower_col + 2, cols - 1)]

    return subset_3x3
