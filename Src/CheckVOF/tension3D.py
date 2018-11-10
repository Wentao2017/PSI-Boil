import sys
noise_file = "uvw-c-press_pall_000858.plt"
 
OpenDatabase(noise_file)

# Contour
AddPlot("Contour", "A")
c = ContourAttributes()
c.colorType = 0
c.contourMethod = 1
c.contourValue = 0.5
c.singleColor = (255, 255, 255, 255)
SetPlotOptions(c)

# View
v = GetView3D()
v.viewNormal = (0, -1, 0)
v.viewUp = (0, 0, 1)
v.perspective = 0
SetView3D(v)

# Vector

DefineVectorExpression("tension", "{U,V,W}")
AddPlot("Vector", "tension")
v = VectorAttributes()
v.useStride = 1
v.stride = 1
SetPlotOptions(v)

DrawPlots()
SaveWindow()
sys.exit() 
