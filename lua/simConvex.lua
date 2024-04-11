local simConvex = loadPlugin('simConvex');

function simConvex.hull(handles)
    local vert = {}
    for _, h in ipairs(handles) do
        local t = sim.getObjectType(h)
        if t == sim.object_shape_type then
            local v = sim.getShapeMesh(h)
            local m = sim.getObjectMatrix(h)
            v = sim.multiplyVector(m, v)
            for _, x in ipairs(v) do table.insert(vert, x) end
        elseif t == sim.object_dummy_type then
            local p = sim.getObjectPosition(h)
            for _, x in ipairs(p) do table.insert(vert, x) end
        elseif t == sim.object_pointcloud_type then
            local v = sim.getPointCloudPoints(h)
            local m = sim.getObjectMatrix(h)
            v = sim.multiplyVector(m, v)
            for _, x in ipairs(v) do table.insert(vert, x) end
        elseif t == sim.object_octree_type then
            local v = sim.getOctreeVoxels(h)
            local vsh = 0.5 * sim.getObjectFloatParam(h, sim.octreefloatparam_voxelsize)
            local m = sim.getObjectMatrix(h)
            local vx = Vector{m[1], m[5], m[9]} * vsh
            local vy = Vector{m[2], m[6], m[10]} * vsh
            local vz = Vector{m[3], m[7], m[11]} * vsh
            for x = -1, 1, 2 do
                for y = -1, 1, 2 do
                    for z = -1, 1, 2 do
                        local d = vx * x + vy * y + vz * z
                        local m2 = table.clone(m)
                        m2[4] = m2[4] + d[1]
                        m2[8] = m2[8] + d[2]
                        m2[12] = m2[12] + d[3]
                        local v2 = sim.multiplyVector(m2, v)
                        for _, xx in ipairs(v2) do table.insert(vert, xx) end
                    end
                end
            end
        else
            sim.addLog(sim.verbosity_warnings, 'unsupported object type: ' .. t)
        end
    end
    if #vert == 0 then error('empty input') end
    local v, i = sim.getQHull(vert)
    local h = sim.createShape(0, 0.0, v, i)
    return h
end

function simConvex.qhull(vertices)
    local vert, ind = sim.getQHull(vertices)
    return vert, ind
end

return simConvex
