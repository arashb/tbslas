from vis_utils import *

def draw_slice():
    AddPlot("Mesh", "mesh", 1, 1)

    AddOperator("Slice", 1)
    SliceAtts = SliceAttributes()
    SliceAtts.originType = SliceAtts.Intercept  # Point, Intercept, Percent, Zone, Node
    SliceAtts.originPoint = (0, 0, 0)
    SliceAtts.originIntercept = 0.5
    SliceAtts.originPercent = 0
    SliceAtts.originZone = 0
    SliceAtts.originNode = 0
    SliceAtts.normal = (0, 0, 1)
    SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
    SliceAtts.upAxis = (0, 1, 0)
    SliceAtts.project2d = 1
    SliceAtts.interactive = 0
    SliceAtts.flip = 0
    SliceAtts.originZoneDomain = 0
    SliceAtts.originNodeDomain = 0
    SliceAtts.meshName = "mesh"
    SliceAtts.theta = 0
    SliceAtts.phi = 90
    SetOperatorOptions(SliceAtts, 1)

    AddPlot("Pseudocolor", "cheb_value", 1, 1)

    DrawPlots()


