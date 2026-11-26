# simConvex API reference

## Alphabetical list

- [simConvex.hacd](#hacd)
- [simConvex.hull](#hull)
- [simConvex.qhull](#qhull)
- [simConvex.vhacd](#vhacd)

## Commands

### <a name="hacd" id="hacd"></a> <a name="simConvex.hacd" id="simConvex.hacd"></a> simConvex.hacd

Computes the HACD convex decomposition of a shape

#### Synopsis:

```lua
int[] convexShapeHandles = simConvex.hacd(int shapeHandle, map params)
```

```python
list convexShapeHandles = simConvex.hacd(int shapeHandle, dict params)
```

#### Parameters:

- **shapeHandle** (int): the handle of the shape
- **params** (map): a map with optional parameters. Can be nil. Allowed parameters are:
  - min_cluster_cnt (default: 1)
  - max_concavity (default: 100.0)
  - max_connection_dist (default: 30.0)
  - triangle_cnt_decimated_mesh (default: 500)
  - max_vertices_cnt (default: 200)
  - small_cluster_detect_threshold (default: 0.25)
  - add_extra_pts (default: true)
  - add_extra_face_pts (default: true)

#### Return values:

- **convexShapeHandles** (int[]): the convex decomposition

#### See also:

[vhacd](#vhacd), [hull](#hull)

### <a name="hull" id="hull"></a> <a name="simConvex.hull" id="simConvex.hull"></a> simConvex.hull

Computes the convex hull of one or more objects

#### Synopsis:

```lua
int convexShapeHandle = simConvex.hull(int[] objectHandles, float growth = 0.0)
```

```python
int convexShapeHandle = simConvex.hull(list objectHandles, float growth = 0.0)
```

#### Parameters:

- **objectHandles** (int[]): the handles of objects. Only shapes, dummies, OC-trees and point clouds are taken into account
- **growth** (float, default: 0.0): an optional growth

#### Return values:

- **convexShapeHandle** (int): the convex hull

#### See also:

[qhull](#qhull), [hacd](#hacd), [vhacd](#vhacd)

### <a name="qhull" id="qhull"></a> <a name="simConvex.qhull" id="simConvex.qhull"></a> simConvex.qhull

Computes the convex hull from 3D points

#### Synopsis:

```lua
float[] vertices, int[] indices = simConvex.qhull(float[] points, float growth = 0.0)
```

```python
list vertices, list indices = simConvex.qhull(list points, float growth = 0.0)
```

#### Parameters:

- **points** (float[]): an array with consecutive 3D point coordinates
- **growth** (float, default: 0.0): an optional growth

#### Return values:

- **vertices** (float[]): the vertices of the convex hull
- **indices** (int[]): the indices of the convex hull

#### See also:

[hull](#hull)

### <a name="vhacd" id="vhacd"></a> <a name="simConvex.vhacd" id="simConvex.vhacd"></a> simConvex.vhacd

Computes the V-HACD convex decomposition of a shape

#### Synopsis:

```lua
int[] convexShapeHandles = simConvex.vhacd(int shapeHandle, map params)
```

```python
list convexShapeHandles = simConvex.vhacd(int shapeHandle, dict params)
```

#### Parameters:

- **shapeHandle** (int): the handle of the shape
- **params** (map): a map with optional parameters. Can be nil. Allowed parameters are:
  - resolution (default: 100000)
  - concavity (default: 0.1)
  - plane_downsampling (default: 4)
  - hull_downsampling (default: 4)
  - alpha (default: 0.05)
  - beta (default: 0.05)
  - max_vertices (default: 64)
  - min_volume (default: 0.0001)
  - pca (default: false)
  - voxels (default: true)

#### Return values:

- **convexShapeHandles** (int[]): the convex decomposition

#### See also:

[hacd](#hacd), [hull](#hull)
