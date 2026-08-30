sub defineGlobalVariable(oProject, varName, varValue)
oProject.ChangeProperty Array("NAME:AllTabs", Array("NAME:ProjectVariableTab", Array("NAME:PropServers",  _
  "ProjectVariables"), Array("NAME:NewProps", Array("NAME:$"+varName, "PropType:=", "VariableProp", "UserDef:=",  _
  true, "Value:=", varValue))))
end sub

sub makeGBHidden(oProject, varName)
oProject.ChangeProperty Array("NAME:AllTabs", Array("NAME:ProjectVariableTab", Array("NAME:PropServers",  _
  "ProjectVariables"), Array("NAME:ChangedProps", Array("NAME:$"+varName, "Hidden:=", true))))
end sub

sub uniteEdges(oEditor, eNames)

oEditor.Unite Array("NAME:Selections", "Selections:=",  _
  eNames), Array("NAME:UniteParameters", "KeepOriginals:=",  _
  false)

end sub 

sub coverLoop(oEditor, lName)
oEditor.CoverLines Array("NAME:Selections", "Selections:=", lName, "NewPartsModelFlag:=",  _
  "Model")
end sub

sub rename(oEditor, oldName, newName)
oEditor.ChangeProperty Array("NAME:AllTabs", Array("NAME:Geometry3DAttributeTab", Array("NAME:PropServers",  _
  oldName), Array("NAME:ChangedProps", Array("NAME:Name", "Value:=", newName))))
end sub

sub subtract(oEditor, toolParts, blankParts)
oEditor.Subtract Array("NAME:Selections", "Blank Parts:=", blankParts, "Tool Parts:=",  _
  toolParts), Array("NAME:SubtractParameters", "KeepOriginals:=", false)
end sub

sub changeObjectColor(oEditor, oName, R, G, B)
oEditor.ChangeProperty Array("NAME:AllTabs", Array("NAME:Geometry3DAttributeTab", Array("NAME:PropServers",  _
  oName), Array("NAME:ChangedProps", Array("NAME:Color", "R:=", R, "G:=", G, "B:=", B))))
end sub

sub drawSegment(oEditor, index1, index2, name)

oEditor.CreatePolyline Array("NAME:PolylineParameters", "IsPolylineCovered:=", true, "IsPolylineClosed:=",  _
  false, Array("NAME:PolylinePoints", Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index1-1)+"]", "Y:=", "$y_pts["+cstr(index1-1)+"]", "Z:=",  _
  "0mm"), Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index2-1)+"]", "Y:=", "$y_pts["+cstr(index2-1)+"]", "Z:=", "0mm")), Array("NAME:PolylineSegments", Array("NAME:PLSegment", "SegmentType:=",  _
  "Line", "StartIndex:=", 0, "NoOfPoints:=", 2)), Array("NAME:PolylineXSection", "XSectionType:=",  _
  "None", "XSectionOrient:=", "Auto", "XSectionWidth:=", "0mm", "XSectionTopWidth:=",  _
  "0mm", "XSectionHeight:=", "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=",  _
  "Corner")), Array("NAME:Attributes", "Name:=", name, "Flags:=", "", "Color:=",  _
  "(143 175 143)", "Transparency:=", 0, "PartCoordinateSystem:=", "Global", "UDMId:=",  _
  "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SurfaceMaterialValue:=",  _
  "" & Chr(34) & "" & Chr(34) & "", "SolveInside:=", true, "ShellElement:=",  _
  false, "ShellElementThickness:=", "0mm", "IsMaterialEditable:=", true, "UseMaterialAppearance:=",  _
  false, "IsLightweight:=", false)

end sub

sub drawArcCPA(oEditor, index1, index2, index3, name)

oEditor.CreatePolyline Array("NAME:PolylineParameters", "IsPolylineCovered:=", true, "IsPolylineClosed:=",  _
  false, Array("NAME:PolylinePoints", Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index2-1)+"]", "Y:=", "$y_pts["+cstr(index2-1)+"]", "Z:=",  _
  "0mm"), Array("NAME:PLPoint", "X:=", "1mm", "Y:=",  _
  "1mm", "Z:=", "0mm"), Array("NAME:PLPoint", "X:=",  _
  "-0.368087051817776mm", "Y:=", "-0.206384113560369mm", "Z:=", "0mm")), Array("NAME:PolylineSegments", Array("NAME:PLSegment", "SegmentType:=",  _
  "AngularArc", "StartIndex:=", 0, "NoOfPoints:=", 3, "NoOfSegments:=", "0", "ArcAngle:=",  _
  "$e_angles["+cstr(index3-1)+"]", "ArcCenterX:=", "$x_pts["+cstr(index1-1)+"]", "ArcCenterY:=", "$y_pts["+cstr(index1-1)+"]", "ArcCenterZ:=",  _
  "0mm", "ArcPlane:=", "XY")), Array("NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=",  _
  "Auto", "XSectionWidth:=", "0mm", "XSectionTopWidth:=", "0mm", "XSectionHeight:=",  _
  "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner")), Array("NAME:Attributes", "Name:=",  _
  name, "Flags:=", "", "Color:=", "(143 175 143)", "Transparency:=", 0, "PartCoordinateSystem:=",  _
  "Global", "UDMId:=", "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SurfaceMaterialValue:=",  _
  "" & Chr(34) & "" & Chr(34) & "", "SolveInside:=", true, "ShellElement:=",  _
  false, "ShellElementThickness:=", "0mm", "IsMaterialEditable:=", true, "UseMaterialAppearance:=",  _
  false, "IsLightweight:=", false)

end sub

sub drawArc(oEditor, index1, index2, index3, name)

arg1 = "atan2(" + "$y_pts["+cstr(index3-1)+"]" + "-" + "$y_pts["+cstr(index1-1)+"]," + "$x_pts["+cstr(index3-1)+"]" + "-" + "$x_pts["+cstr(index1-1)+"])"
arg2 = "atan2(" + "$y_pts["+cstr(index2-1)+"]" + "-" + "$y_pts["+cstr(index1-1)+"]," + "$x_pts["+cstr(index2-1)+"]" + "-" + "$x_pts["+cstr(index1-1)+"])"

oEditor.CreatePolyline Array("NAME:PolylineParameters", "IsPolylineCovered:=", true, "IsPolylineClosed:=",  _
  false, Array("NAME:PolylinePoints", Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index2-1)+"]", "Y:=", "$y_pts["+cstr(index2-1)+"]", "Z:=",  _
  "0mm"), Array("NAME:PLPoint", "X:=", "1mm", "Y:=",  _
  "1mm", "Z:=", "0mm"), Array("NAME:PLPoint", "X:=",  _
  "-0.368087051817776mm", "Y:=", "-0.206384113560369mm", "Z:=", "0mm")), Array("NAME:PolylineSegments", Array("NAME:PLSegment", "SegmentType:=",  _
  "AngularArc", "StartIndex:=", 0, "NoOfPoints:=", 3, "NoOfSegments:=", "0", "ArcAngle:=",  _
  arg1+"-"+arg2, "ArcCenterX:=", "$x_pts["+cstr(index1-1)+"]", "ArcCenterY:=", "$y_pts["+cstr(index1-1)+"]", "ArcCenterZ:=",  _
  "0mm", "ArcPlane:=", "XY")), Array("NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=",  _
  "Auto", "XSectionWidth:=", "0mm", "XSectionTopWidth:=", "0mm", "XSectionHeight:=",  _
  "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner")), Array("NAME:Attributes", "Name:=",  _
  name, "Flags:=", "", "Color:=", "(143 175 143)", "Transparency:=", 0, "PartCoordinateSystem:=",  _
  "Global", "UDMId:=", "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SurfaceMaterialValue:=",  _
  "" & Chr(34) & "" & Chr(34) & "", "SolveInside:=", true, "ShellElement:=",  _
  false, "ShellElementThickness:=", "0mm", "IsMaterialEditable:=", true, "UseMaterialAppearance:=",  _
  false, "IsLightweight:=", false)

end sub

Dim oAnsoftApp
Dim oDesktop
Dim oProject
Dim oDesign
Dim oEditor
Dim oModule
Set oAnsoftApp = CreateObject("Ansoft.ElectronicsDesktop")
Set oDesktop = oAnsoftApp.GetAppDesktop()
oDesktop.RestoreWindow
Set oProject = oDesktop.NewProject
oProject.InsertDesign "Maxwell 2D", "NewDesign", "Magnetostatic", ""
Set oDesign = oProject.SetActiveDesign("NewDesign")
Set oEditor = oDesign.SetActiveEditor("3D Modeler")

call defineGlobalVariable(oProject, "x_pts", "[0,44.9937495659,45.9936106674,46.313212429,56.7706554351,60,72.5,72.2241156117,44.8287614141,45,46.3583063575,57.0295474348,20,27.0479012879,30.512002903,44.3,38.3649253877,17.3205080757,36.0718519318,39.5359535469,37.8039027394,27.0479012879,36.677208179,39.5303777279,42.2694910946,39.4163215457,36.677208179,40.8429063201,19.918584287,22.9496732003,25.9807621135,24.6996732785,38.8039027394,41.3064318117] mm")
call makeGBHidden(oProject, "x_pts")
call defineGlobalVariable(oProject, "y_pts", "[0,0.75,0.766666666667,2.04424138392,2.95914909642,0,0,6.31879134921,3.92200842364,0,0,0,0,2,0,0,22.15,10,17.6299410003,15.6299410003,16.6299410003,0,0,0,8.43012411501,9.35717509813,0.927050983125,8.89364960657,11.5,13.25,15,10.2189111319,18.3619918079,10.320234381] mm")
call makeGBHidden(oProject, "y_pts")
call defineGlobalVariable(oProject, "e_angles", "[0,0,0,-95,0,5,0,-4.04502612622,0,0,0.954973873785,2.52736574935,0,0,0,0,0,30,0,-30,0,0,0,0,90,0,0,0,0,0,90,0,0,0,90.0000014788,89.9999985212,90,90] deg")
call makeGBHidden(oProject, "e_angles")
call drawSegment(oEditor, 2, 3, "stator_loop_1_1")
call drawSegment(oEditor, 3, 4, "stator_loop_1_2")
call drawSegment(oEditor, 4, 5, "stator_loop_1_3")
call drawArcCPA(oEditor, 12, 5, 4, "stator_loop_1_4")
call drawSegment(oEditor, 6, 7, "stator_loop_1_5")
call drawArcCPA(oEditor, 1, 7, 6, "stator_loop_1_6")
call drawSegment(oEditor, 8, 9, "stator_loop_1_7")
call drawArcCPA(oEditor, 1, 9, 8, "stator_loop_1_8")
call uniteEdges(oEditor, "stator_loop_1_1,stator_loop_1_2,stator_loop_1_3,stator_loop_1_4,stator_loop_1_5,stator_loop_1_6,stator_loop_1_7,stator_loop_1_8")
call coverLoop(oEditor, "stator_loop_1_1")
call rename(oEditor, "stator_loop_1_1", "stator")
call changeObjectColor(oEditor, "stator", 200,200,200)
call drawSegment(oEditor, 11, 6, "sc_loop_1_10")
call drawArcCPA(oEditor, 12, 5, 4, "sc_loop_1_4")
call drawSegment(oEditor, 4, 5, "sc_loop_1_3")
call drawArcCPA(oEditor, 1, 11, 12, "sc_loop_1_12")
call uniteEdges(oEditor, "sc_loop_1_10,sc_loop_1_4,sc_loop_1_3,sc_loop_1_12")
call coverLoop(oEditor, "sc_loop_1_10")
call rename(oEditor, "sc_loop_1_10", "sc")
call changeObjectColor(oEditor, "sc", 255,137,39)
call drawSegment(oEditor, 10, 11, "sap_loop_1_9")
call drawArcCPA(oEditor, 1, 11, 12, "sap_loop_1_12")
call drawSegment(oEditor, 3, 4, "sap_loop_1_2")
call drawSegment(oEditor, 2, 3, "sap_loop_1_1")
call drawArcCPA(oEditor, 1, 10, 11, "sap_loop_1_11")
call uniteEdges(oEditor, "sap_loop_1_9,sap_loop_1_12,sap_loop_1_2,sap_loop_1_1,sap_loop_1_11")
call coverLoop(oEditor, "sap_loop_1_9")
call rename(oEditor, "sap_loop_1_9", "sap")
call changeObjectColor(oEditor, "sap", 0,255,255)
call drawSegment(oEditor, 13, 22, "rotor_loop_1_13")
call drawSegment(oEditor, 14, 22, "rotor_loop_1_26")
call drawSegment(oEditor, 19, 14, "rotor_loop_1_24")
call drawArcCPA(oEditor, 21, 33, 37, "rotor_loop_1_37")
call drawArcCPA(oEditor, 21, 20, 25, "rotor_loop_1_25")
call drawSegment(oEditor, 15, 20, "rotor_loop_1_22")
call drawSegment(oEditor, 15, 23, "rotor_loop_1_15")
call drawSegment(oEditor, 27, 23, "rotor_loop_1_32")
call drawSegment(oEditor, 26, 27, "rotor_loop_1_29")
call drawArcCPA(oEditor, 28, 34, 38, "rotor_loop_1_38")
call drawArcCPA(oEditor, 28, 25, 31, "rotor_loop_1_31")
call drawSegment(oEditor, 24, 25, "rotor_loop_1_27")
call drawSegment(oEditor, 24, 16, "rotor_loop_1_17")
call drawArcCPA(oEditor, 1, 16, 18, "rotor_loop_1_18")
call drawSegment(oEditor, 17, 31, "rotor_loop_1_19")
call drawArcCPA(oEditor, 30, 32, 36, "rotor_loop_1_36")
call drawArcCPA(oEditor, 30, 29, 35, "rotor_loop_1_35")
call drawSegment(oEditor, 29, 18, "rotor_loop_1_34")
call drawArcCPA(oEditor, 1, 18, 20, "rotor_loop_1_20")
call uniteEdges(oEditor, "rotor_loop_1_13,rotor_loop_1_26,rotor_loop_1_24,rotor_loop_1_37,rotor_loop_1_25,rotor_loop_1_22,rotor_loop_1_15,rotor_loop_1_32,rotor_loop_1_29,rotor_loop_1_38,rotor_loop_1_31,rotor_loop_1_27,rotor_loop_1_17,rotor_loop_1_18,rotor_loop_1_19,rotor_loop_1_36,rotor_loop_1_35,rotor_loop_1_34,rotor_loop_1_20")
call coverLoop(oEditor, "rotor_loop_1_13")
call rename(oEditor, "rotor_loop_1_13", "rotor")
call changeObjectColor(oEditor, "rotor", 200,200,200)
call drawSegment(oEditor, 14, 15, "magnet1_loop_1_21")
call drawSegment(oEditor, 15, 20, "magnet1_loop_1_22")
call drawSegment(oEditor, 20, 19, "magnet1_loop_1_23")
call drawSegment(oEditor, 19, 14, "magnet1_loop_1_24")
call uniteEdges(oEditor, "magnet1_loop_1_21,magnet1_loop_1_22,magnet1_loop_1_23,magnet1_loop_1_24")
call coverLoop(oEditor, "magnet1_loop_1_21")
call rename(oEditor, "magnet1_loop_1_21", "magnet1")
call changeObjectColor(oEditor, "magnet1", 160,78,146)
call drawSegment(oEditor, 24, 25, "magnet2_loop_1_27")
call drawSegment(oEditor, 25, 26, "magnet2_loop_1_28")
call drawSegment(oEditor, 26, 27, "magnet2_loop_1_29")
call drawSegment(oEditor, 27, 24, "magnet2_loop_1_30")
call uniteEdges(oEditor, "magnet2_loop_1_27,magnet2_loop_1_28,magnet2_loop_1_29,magnet2_loop_1_30")
call coverLoop(oEditor, "magnet2_loop_1_27")
call rename(oEditor, "magnet2_loop_1_27", "magnet2")
call changeObjectColor(oEditor, "magnet2", 160,78,146)
call drawSegment(oEditor, 20, 19, "rap1_loop_1_23")
call drawArcCPA(oEditor, 21, 20, 25, "rap1_loop_1_25")
call drawArcCPA(oEditor, 21, 33, 37, "rap1_loop_1_37")
call uniteEdges(oEditor, "rap1_loop_1_23,rap1_loop_1_25,rap1_loop_1_37")
call coverLoop(oEditor, "rap1_loop_1_23")
call rename(oEditor, "rap1_loop_1_23", "rap1")
call changeObjectColor(oEditor, "rap1", 0,255,255)
call drawSegment(oEditor, 25, 26, "rap2_loop_1_28")
call drawArcCPA(oEditor, 28, 25, 31, "rap2_loop_1_31")
call drawArcCPA(oEditor, 28, 34, 38, "rap2_loop_1_38")
call uniteEdges(oEditor, "rap2_loop_1_28,rap2_loop_1_31,rap2_loop_1_38")
call coverLoop(oEditor, "rap2_loop_1_28")
call rename(oEditor, "rap2_loop_1_28", "rap2")
call changeObjectColor(oEditor, "rap2", 0,255,255)
call drawSegment(oEditor, 22, 15, "rap3_loop_1_14")
call drawSegment(oEditor, 14, 15, "rap3_loop_1_21")
call drawSegment(oEditor, 14, 22, "rap3_loop_1_26")
call uniteEdges(oEditor, "rap3_loop_1_14,rap3_loop_1_21,rap3_loop_1_26")
call coverLoop(oEditor, "rap3_loop_1_14")
call rename(oEditor, "rap3_loop_1_14", "rap3")
call changeObjectColor(oEditor, "rap3", 0,255,255)
call drawSegment(oEditor, 23, 24, "rap4_loop_1_16")
call drawSegment(oEditor, 27, 24, "rap4_loop_1_30")
call drawSegment(oEditor, 27, 23, "rap4_loop_1_32")
call uniteEdges(oEditor, "rap4_loop_1_16,rap4_loop_1_30,rap4_loop_1_32")
call coverLoop(oEditor, "rap4_loop_1_16")
call rename(oEditor, "rap4_loop_1_16", "rap4")
call changeObjectColor(oEditor, "rap4", 0,255,255)
call drawSegment(oEditor, 31, 29, "rap5_loop_1_33")
call drawArcCPA(oEditor, 30, 29, 35, "rap5_loop_1_35")
call drawArcCPA(oEditor, 30, 32, 36, "rap5_loop_1_36")
call uniteEdges(oEditor, "rap5_loop_1_33,rap5_loop_1_35,rap5_loop_1_36")
call coverLoop(oEditor, "rap5_loop_1_33")
call rename(oEditor, "rap5_loop_1_33", "rap5")
call changeObjectColor(oEditor, "rap5", 0,255,255)

