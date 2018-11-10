import sys
noise_file = "uvw-c-press_pall_000004.plt"
 
OpenDatabase(noise_file)

# Contour
AddPlot("Contour", "A")
AddOperator("Slice")
c = ContourAttributes()
c.colorType = 0
c.contourMethod = 1
c.contourValue = 0.5
c.singleColor = (0, 0, 0, 255)
SetPlotOptions(c)

# Vector
DefineVectorExpression("tension", "{U,V,W}")
AddPlot("Vector", "tension")
v = VectorAttributes()
v.useStride = 1
v.stride = 2
SetPlotOptions(v)
AddOperator("Slice")

DrawPlots()
SaveWindow()
sys.exit() 
