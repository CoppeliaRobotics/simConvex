local codeEditorInfos = [[
int[] convexShapeHandles = simConvex.hacd(int shapeHandle, map params)
int[] convexShapeHandles = simConvex.vhacd(int shapeHandle, map params)
int convexShapeHandle = simConvex.hull(int[] objectHandles)
float[] vertices, int[] indices = simConvex.qhull(float[] points)
]]

registerCodeEditorInfos("simConvex", codeEditorInfos)
