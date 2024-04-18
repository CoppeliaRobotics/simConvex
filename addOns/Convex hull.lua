function sysCall_info()
    return {autoStart = false, menu = 'Geometry / Mesh\nConvex hull...'}
end

function sysCall_nonSimulation()
    if leaveNow then
        simUI.destroy(ui)
        if params then
            local generated = {}
            if params.mode == 'perObject' or params.mode == 'perComponent' then
                for i = 1, #params.sel do
                    local obj = params.sel[i]
                    local t = sim.getObjectType(obj)
                    if t == sim.object_shape_type or t == sim.object_pointcloud_type or t == sim.object_octree_type then
                        if t == sim.object_shape_type and params.mode == 'perComponent' then
                            local tt = sim.getShapeGeomInfo(obj)
                            if (tt & 1) > 0 then
                                local shapes = extractSimpleShapes(sim.copyPasteObjects({obj}, 2|4|8|16|32))                            
                                for j = 1, #shapes do
                                    local h = simConvex.hull({shapes[j]}, params.growth)
                                    generated[#generated + 1] = h
                                end
                                sim.removeObjects(shapes)
                            else
                                local h = simConvex.hull({obj}, params.growth)
                                generated[#generated + 1] = h
                            end
                        else
                            local h = simConvex.hull({obj}, params.growth)
                            generated[#generated + 1] = h
                        end
                    end
                end
            else
                -- perSelection
                if #params.sel > 0 then
                    local h = simConvex.hull(params.sel, params.growth)
                    generated[#generated + 1] = h
                end
            end
            local convex = true
            local faces = 0
            for i = 1, #generated do
                local h = generated[i]
                if (sim.getShapeGeomInfo(h) & 4) == 0 then
                    convex = false
                end
                -- Pose, BB:
                sim.relocateShapeFrame(h, {0, 0, 0, 0, 0, 0, 0})
                sim.alignShapeBB(h, {0, 0, 0, 0, 0, 0, 0})

                -- Various:
                sim.setObjectAlias(h, 'convexHull')
                sim.setObjectFloatParam(h, sim.shapefloatparam_shading_angle, 45.0 * math.pi / 180.0)
                local vert, ind = sim.getShapeMesh(h)
                faces = faces + #ind / 3
            end
            if not convex then
                sim.addLog(sim.verbosity_scripterrors, 'One or more of the generated shapes is not convex.') 
            end
            sim.setObjectSel(generated)

            if #generated > 0 then
                sim.announceSceneContentChange()
                sim.addLog(sim.verbosity_scriptinfos, string.format('Generated %i convex hull(s) (with a total of %i triangular faces)', #generated, faces))
            else
                simUI.msgBox(simUI.msgbox_type.info, simUI.msgbox_buttons.ok, "Convex Hull Generator", 'The resulting selection is effectively empty...')
                sim.setObjectSel({})
            end
        end
        return {cmd = 'cleanup'} 
    end
end
    
function sysCall_init()
    sim = require('sim')
    simUI = require('simUI')
    simConvex = require('simConvex')

    local sel = sim.getObjectSel()
    if #sel == 0 or sim.getSimulationState() ~= sim.simulation_stopped then
        simUI.msgBox(simUI.msgbox_type.info, simUI.msgbox_buttons.ok, "Convex Hull Generator", 'Make sure that at least one object is selected, and that simulation is not running.')
    else
        ui = simUI.create(
          [[<ui title="Convex Hull Generator" closeable="true" on-close="onClose" modal="true">
            <radiobutton text="generate one single hull" on-click="updateUi" id="${perSelection}" checked="true"/>
            <radiobutton text="generate one hull per object" on-click="updateUi" id="${perObject}" checked="false" />
            <radiobutton text="generate one hull per object component" on-click="updateUi" id="${perComponent}" checked="false"/>
            <checkbox id="${model_objects}" text="include model objects" checked="false" on-change="updateUi" />
            <checkbox id="${hidden_objects}" text="exclude hidden objects" checked="false" on-change="updateUi" />
            <group flat="true" content-margins="0,0,0,0" layout="form">
                <label text="growth:" />
                <spinbox id="${growth}" minimum="0.0" maximum="100.0" value="0.0" step="0.02" on-change="updateUi" />
            </group>
            <button id="${gen}" text="Generate" on-click="initGenerate" />
        </ui>]]
             )
    end
end

function onClose()
    leaveNow = true
    abort = true
end

function updateUi()
end

function initGenerate()
    local includeModelObjects = simUI.getCheckboxValue(ui, model_objects) > 0
    local excludeHiddenObjects = simUI.getCheckboxValue(ui, hidden_objects) > 0
    local s = sim.getObjectSel()
    local selMap = {}
    for i = 1, #s do
        local h = s[i]
        if sim.getModelProperty(h) == sim.modelproperty_not_model or not includeModelObjects then
            selMap[h] = true
        else
            local tree = sim.getObjectsInTree(h, sim.object_shape_type)
            for j = 1, #tree do
                selMap[tree[j]] = true
            end
        end
    end
    local sel = {}
    for obj, v in pairs(selMap) do
        if not excludeHiddenObjects or (sim.getObjectInt32Param(obj, sim.objintparam_visible) > 0) then
            sel[#sel + 1] = obj
        end
    end
    
    leaveNow = true
    params = {sel = sel}
    
    if simUI.getRadiobuttonValue(ui, perSelection) == 1 then
        params.mode = 'perSelection'
    end
    if simUI.getRadiobuttonValue(ui, perObject) == 1 then
        params.mode = 'perObject'
    end
    if simUI.getRadiobuttonValue(ui, perComponent) == 1 then
        params.mode = 'perComponent'
    end
    params.growth = tonumber(simUI.getSpinboxValue(ui, growth))
end

function extractSimpleShapes(shapes)
    local retVal = {}
    for i = 1, #shapes do
        local shape = shapes[i]
        local t = sim.getShapeGeomInfo(shape)
        if t & 1 > 0 then
            local nshapes = sim.ungroupShape(shape)
            retVal = table.add(retVal, extractSimpleShapes(nshapes))
        else
            retVal[#retVal + 1] = shape
        end
    end
    return retVal
end
