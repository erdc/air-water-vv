try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

quiescent_water_test_gauges_p0_xmf = XDMFReader( FileName='\\\\Hrw-uk.local\\projects\\live\\cas1215$\\model\\Proteus Validation\\quiescent_water_test_gauges\\quiescent_water_test_gauges_p0.xmf' )

AnimationScene1 = GetAnimationScene()
quiescent_water_test_gauges_p0_xmf.CellArrays = ['CellMapL2G', 'elementMaterialTypes']
quiescent_water_test_gauges_p0_xmf.Sets = []
quiescent_water_test_gauges_p0_xmf.PointArrays = ['NodeMapL2G', 'nodeMaterialTypes', 'p', 'phi', 'phiCorr', 'phid', 'u', 'v', 'velocity', 'vof']

AnimationScene1.PlayMode = 'Snap To TimeSteps'

quiescent_water_test_gauges_p1_xmf = XDMFReader( FileName='\\\\Hrw-uk.local\\projects\\live\\cas1215$\\model\\Proteus Validation\\quiescent_water_test_gauges\\quiescent_water_test_gauges_p1.xmf' )

quiescent_water_test_gauges_p1_xmf.CellArrays = ['CellMapL2G', 'elementMaterialTypes']
quiescent_water_test_gauges_p1_xmf.Sets = []
quiescent_water_test_gauges_p1_xmf.PointArrays = ['NodeMapL2G', 'nodeMaterialTypes', 'p', 'phi', 'phiCorr', 'phid', 'u', 'v', 'velocity', 'vof']

quiescent_water_test_gauges_p0_xmf.Grids = ['Grid_614', 'Grid_617', 'Grid_620', 'Grid_623', 'Grid_626', 'Grid_629', 'Grid_632', 'Grid_635', 'Grid_638', 'Grid_641', 'Grid_644', 'Grid_647', 'Grid_650', 'Grid_653', 'Grid_656', 'Grid_659', 'Grid_662', 'Grid_665', 'Grid_668', 'Grid_671', 'Grid_674', 'Grid_677', 'Grid_680', 'Grid_683', 'Grid_686', 'Grid_689', 'Grid_692', 'Grid_695', 'Grid_698', 'Grid_701', 'Grid_704', 'Grid_707', 'Grid_710', 'Grid_713', 'Grid_716', 'Grid_719', 'Grid_722', 'Grid_725', 'Grid_728', 'Grid_731', 'Grid_734', 'Grid_737', 'Grid_740', 'Grid_743', 'Grid_746', 'Grid_749', 'Grid_752', 'Grid_755', 'Grid_758', 'Grid_761', 'Grid_764', 'Grid_767', 'Grid_770', 'Grid_773', 'Grid_776', 'Grid_779', 'Grid_782', 'Grid_785', 'Grid_788', 'Grid_791', 'Grid_794', 'Grid_797', 'Grid_800', 'Grid_803', 'Grid_806', 'Grid_809', 'Grid_812', 'Grid_815', 'Grid_818', 'Grid_821', 'Grid_824', 'Grid_827', 'Grid_830', 'Grid_833', 'Grid_836', 'Grid_839', 'Grid_842', 'Grid_845', 'Grid_848', 'Grid_851', 'Grid_854', 'Grid_857', 'Grid_860', 'Grid_863', 'Grid_866', 'Grid_869', 'Grid_872', 'Grid_875', 'Grid_878', 'Grid_881', 'Grid_884', 'Grid_887', 'Grid_890', 'Grid_893', 'Grid_896', 'Grid_899', 'Grid_902', 'Grid_905', 'Grid_908', 'Grid_911', 'Grid_914', 'Grid_917']

quiescent_water_test_gauges_p1_xmf.Grids = ['Grid_920', 'Grid_923', 'Grid_926', 'Grid_929', 'Grid_932', 'Grid_935', 'Grid_938', 'Grid_941', 'Grid_944', 'Grid_947', 'Grid_950', 'Grid_953', 'Grid_956', 'Grid_959', 'Grid_962', 'Grid_965', 'Grid_968', 'Grid_971', 'Grid_974', 'Grid_977', 'Grid_980', 'Grid_983', 'Grid_986', 'Grid_989', 'Grid_992', 'Grid_995', 'Grid_998', 'Grid_1001', 'Grid_1004', 'Grid_1007', 'Grid_1010', 'Grid_1013', 'Grid_1016', 'Grid_1019', 'Grid_1022', 'Grid_1025', 'Grid_1028', 'Grid_1031', 'Grid_1034', 'Grid_1037', 'Grid_1040', 'Grid_1043', 'Grid_1046', 'Grid_1049', 'Grid_1052', 'Grid_1055', 'Grid_1058', 'Grid_1061', 'Grid_1064', 'Grid_1067', 'Grid_1070', 'Grid_1073', 'Grid_1076', 'Grid_1079', 'Grid_1082', 'Grid_1085', 'Grid_1088', 'Grid_1091', 'Grid_1094', 'Grid_1097', 'Grid_1100', 'Grid_1103', 'Grid_1106', 'Grid_1109', 'Grid_1112', 'Grid_1115', 'Grid_1118', 'Grid_1121', 'Grid_1124', 'Grid_1127', 'Grid_1130', 'Grid_1133', 'Grid_1136', 'Grid_1139', 'Grid_1142', 'Grid_1145', 'Grid_1148', 'Grid_1151', 'Grid_1154', 'Grid_1157', 'Grid_1160', 'Grid_1163', 'Grid_1166', 'Grid_1169', 'Grid_1172', 'Grid_1175', 'Grid_1178', 'Grid_1181', 'Grid_1184', 'Grid_1187', 'Grid_1190', 'Grid_1193', 'Grid_1196', 'Grid_1199', 'Grid_1202', 'Grid_1205', 'Grid_1208', 'Grid_1211', 'Grid_1214', 'Grid_1217', 'Grid_1220', 'Grid_1223']

RenderView1 = GetRenderView()
a1_NodeMapL2G_PVLookupTable = GetLookupTableForArray( "NodeMapL2G", 1 )

SetActiveSource(quiescent_water_test_gauges_p0_xmf)
DataRepresentation5 = Show()
DataRepresentation5.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation5.SelectionPointFieldDataArrayName = 'NodeMapL2G'
DataRepresentation5.SelectionCellFieldDataArrayName = 'CellMapL2G'
DataRepresentation5.ScalarOpacityFunction = []
DataRepresentation5.ColorArrayName = 'NodeMapL2G'
DataRepresentation5.ScalarOpacityUnitDistance = 0.0970043273963236
DataRepresentation5.LookupTable = a1_NodeMapL2G_PVLookupTable
DataRepresentation5.ScaleFactor = 0.20874056097392868

RenderView1.CenterOfRotation = [1.0437028048696433, 0.9, 0.0]

SetActiveSource(quiescent_water_test_gauges_p1_xmf)
DataRepresentation6 = Show()
DataRepresentation6.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation6.SelectionPointFieldDataArrayName = 'NodeMapL2G'
DataRepresentation6.SelectionCellFieldDataArrayName = 'CellMapL2G'
DataRepresentation6.ScalarOpacityUnitDistance = 0.0965156371233247
DataRepresentation6.ScaleFactor = 0.20174900000000007

RenderView1.CameraPosition = [1.0437028048696433, 0.9, 5.324788113397622]
RenderView1.CameraFocalPoint = [1.0437028048696433, 0.9, 0.0]
RenderView1.CameraClippingRange = [5.2715402322636455, 5.404659935098586]
RenderView1.CameraParallelScale = 1.3781565748828255

DataRepresentation6.ScalarOpacityFunction = []
DataRepresentation6.ColorArrayName = 'NodeMapL2G'
DataRepresentation6.LookupTable = a1_NodeMapL2G_PVLookupTable

GroupDatasets2 = GroupDatasets( Input=[ quiescent_water_test_gauges_p0_xmf, quiescent_water_test_gauges_p1_xmf ] )

DataRepresentation7 = Show()
DataRepresentation7.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation7.SelectionPointFieldDataArrayName = 'NodeMapL2G'
DataRepresentation7.SelectionCellFieldDataArrayName = 'CellMapL2G'
DataRepresentation7.ColorArrayName = 'NodeMapL2G'
DataRepresentation7.ScalarOpacityUnitDistance = 0.10377056257834123
DataRepresentation7.ExtractedBlockIndex = 1
DataRepresentation7.ScaleFactor = 0.32200000000000006

PlotSelectionOverTime2 = PlotSelectionOverTime()

DataRepresentation5.Visibility = 0

DataRepresentation6.Visibility = 0

DataRepresentation7.ScalarOpacityFunction = []
DataRepresentation7.LookupTable = a1_NodeMapL2G_PVLookupTable

selection_source_5367 = LocationSelectionSource( ContainingCells=0, InsideOut=0, FieldType='POINT', Locations=[0.5, 0.5, 0.0] )

PlotSelectionOverTime2.Selection = selection_source_5367

XYChartView2 = CreateXYPlotView()
XYChartView2.ViewTime = 0.0

DataRepresentation8 = Show()
DataRepresentation8.XArrayName = 'Time'
DataRepresentation8.SeriesVisibility = ['velocity (0)', '0', 'velocity (1)', '0', 'velocity (2)', '0', 'vtkValidPointMask', '0', 'Time', '0', 'Probe Coordinates (0)', '0', 'Probe Coordinates (1)', '0', 'Probe Coordinates (2)', '0', 'vtkOriginalIndices', '0']
DataRepresentation8.CompositeDataSetIndex = 1
DataRepresentation8.AttributeType = 'Row Data'
DataRepresentation8.UseIndexForXAxis = 0

AnimationScene1.ViewModules = [ RenderView1, XYChartView2 ]

Render()
