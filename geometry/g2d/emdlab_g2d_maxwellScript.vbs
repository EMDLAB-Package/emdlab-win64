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

call defineGlobalVariable(oProject, "x_pts", "[-246,-212,-212,-202,-202,-194,-194,-184,-184,-176,-176,-166,-166,-158,-158,-148,-148,-140,-140,-130,-130,-122,-122,-112,-112,-104,-104,-94,-94,-86,-86,-76,-76,-68,-68,-58,-58,-50,-50,-40,-40,-32,-32,-22,-22,-14,-14,-4,-4,4,4,14,14,22,22,32,32,40,40,50,50,58,58,68,68,76,76,86,86,94,94,104,104,112,112,122,122,130,130,140,140,148,148,158,158,166,166,176,176,184,184,194,194,202,202,212,212,246,261.541053161,-261.541053161,432,432,-432,-432,-432,432,-432,432,-432,432,-432,432] mm")
call makeGBHidden(oProject, "x_pts")
call defineGlobalVariable(oProject, "y_pts", "[5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,45,45,5,5,63,63,5,98,98,5,-47,-47,-17,-17,-9,-9,-5,-5] mm")
call makeGBHidden(oProject, "y_pts")
call defineGlobalVariable(oProject, "e_angles", "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] deg")
call makeGBHidden(oProject, "e_angles")
call drawSegment(oEditor, 2, 5, "ca_1_loop_1_e101")
call drawSegment(oEditor, 4, 5, "ca_1_loop_1_e4")
call drawSegment(oEditor, 3, 4, "ca_1_loop_1_e3")
call drawSegment(oEditor, 2, 3, "ca_1_loop_1_e2")
call uniteEdges(oEditor, "ca_1_loop_1_e101,ca_1_loop_1_e4,ca_1_loop_1_e3,ca_1_loop_1_e2")
call coverLoop(oEditor, "ca_1_loop_1_e101")
call rename(oEditor, "ca_1_loop_1_e101", "ca_1")
call changeObjectColor(oEditor, "ca_1", 255,127,39)
call drawSegment(oEditor, 6, 9, "ca_2_loop_1_e102")
call drawSegment(oEditor, 8, 9, "ca_2_loop_1_e8")
call drawSegment(oEditor, 7, 8, "ca_2_loop_1_e7")
call drawSegment(oEditor, 6, 7, "ca_2_loop_1_e6")
call uniteEdges(oEditor, "ca_2_loop_1_e102,ca_2_loop_1_e8,ca_2_loop_1_e7,ca_2_loop_1_e6")
call coverLoop(oEditor, "ca_2_loop_1_e102")
call rename(oEditor, "ca_2_loop_1_e102", "ca_2")
call changeObjectColor(oEditor, "ca_2", 255,127,39)
call drawSegment(oEditor, 10, 13, "ca_3_loop_1_e103")
call drawSegment(oEditor, 12, 13, "ca_3_loop_1_e12")
call drawSegment(oEditor, 11, 12, "ca_3_loop_1_e11")
call drawSegment(oEditor, 10, 11, "ca_3_loop_1_e10")
call uniteEdges(oEditor, "ca_3_loop_1_e103,ca_3_loop_1_e12,ca_3_loop_1_e11,ca_3_loop_1_e10")
call coverLoop(oEditor, "ca_3_loop_1_e103")
call rename(oEditor, "ca_3_loop_1_e103", "ca_3")
call changeObjectColor(oEditor, "ca_3", 255,127,39)
call drawSegment(oEditor, 14, 17, "ca_4_loop_1_e104")
call drawSegment(oEditor, 16, 17, "ca_4_loop_1_e16")
call drawSegment(oEditor, 15, 16, "ca_4_loop_1_e15")
call drawSegment(oEditor, 14, 15, "ca_4_loop_1_e14")
call uniteEdges(oEditor, "ca_4_loop_1_e104,ca_4_loop_1_e16,ca_4_loop_1_e15,ca_4_loop_1_e14")
call coverLoop(oEditor, "ca_4_loop_1_e104")
call rename(oEditor, "ca_4_loop_1_e104", "ca_4")
call changeObjectColor(oEditor, "ca_4", 255,127,39)
call drawSegment(oEditor, 18, 21, "ca_5_loop_1_e105")
call drawSegment(oEditor, 20, 21, "ca_5_loop_1_e20")
call drawSegment(oEditor, 19, 20, "ca_5_loop_1_e19")
call drawSegment(oEditor, 18, 19, "ca_5_loop_1_e18")
call uniteEdges(oEditor, "ca_5_loop_1_e105,ca_5_loop_1_e20,ca_5_loop_1_e19,ca_5_loop_1_e18")
call coverLoop(oEditor, "ca_5_loop_1_e105")
call rename(oEditor, "ca_5_loop_1_e105", "ca_5")
call changeObjectColor(oEditor, "ca_5", 255,127,39)
call drawSegment(oEditor, 22, 25, "ca_6_loop_1_e106")
call drawSegment(oEditor, 24, 25, "ca_6_loop_1_e24")
call drawSegment(oEditor, 23, 24, "ca_6_loop_1_e23")
call drawSegment(oEditor, 22, 23, "ca_6_loop_1_e22")
call uniteEdges(oEditor, "ca_6_loop_1_e106,ca_6_loop_1_e24,ca_6_loop_1_e23,ca_6_loop_1_e22")
call coverLoop(oEditor, "ca_6_loop_1_e106")
call rename(oEditor, "ca_6_loop_1_e106", "ca_6")
call changeObjectColor(oEditor, "ca_6", 255,127,39)
call drawSegment(oEditor, 26, 29, "ca_7_loop_1_e107")
call drawSegment(oEditor, 28, 29, "ca_7_loop_1_e28")
call drawSegment(oEditor, 27, 28, "ca_7_loop_1_e27")
call drawSegment(oEditor, 26, 27, "ca_7_loop_1_e26")
call uniteEdges(oEditor, "ca_7_loop_1_e107,ca_7_loop_1_e28,ca_7_loop_1_e27,ca_7_loop_1_e26")
call coverLoop(oEditor, "ca_7_loop_1_e107")
call rename(oEditor, "ca_7_loop_1_e107", "ca_7")
call changeObjectColor(oEditor, "ca_7", 255,127,39)
call drawSegment(oEditor, 30, 33, "ca_8_loop_1_e108")
call drawSegment(oEditor, 32, 33, "ca_8_loop_1_e32")
call drawSegment(oEditor, 31, 32, "ca_8_loop_1_e31")
call drawSegment(oEditor, 30, 31, "ca_8_loop_1_e30")
call uniteEdges(oEditor, "ca_8_loop_1_e108,ca_8_loop_1_e32,ca_8_loop_1_e31,ca_8_loop_1_e30")
call coverLoop(oEditor, "ca_8_loop_1_e108")
call rename(oEditor, "ca_8_loop_1_e108", "ca_8")
call changeObjectColor(oEditor, "ca_8", 255,127,39)
call drawSegment(oEditor, 34, 37, "ca_9_loop_1_e109")
call drawSegment(oEditor, 36, 37, "ca_9_loop_1_e36")
call drawSegment(oEditor, 35, 36, "ca_9_loop_1_e35")
call drawSegment(oEditor, 34, 35, "ca_9_loop_1_e34")
call uniteEdges(oEditor, "ca_9_loop_1_e109,ca_9_loop_1_e36,ca_9_loop_1_e35,ca_9_loop_1_e34")
call coverLoop(oEditor, "ca_9_loop_1_e109")
call rename(oEditor, "ca_9_loop_1_e109", "ca_9")
call changeObjectColor(oEditor, "ca_9", 255,127,39)
call drawSegment(oEditor, 38, 41, "ca_10_loop_1_e110")
call drawSegment(oEditor, 40, 41, "ca_10_loop_1_e40")
call drawSegment(oEditor, 39, 40, "ca_10_loop_1_e39")
call drawSegment(oEditor, 38, 39, "ca_10_loop_1_e38")
call uniteEdges(oEditor, "ca_10_loop_1_e110,ca_10_loop_1_e40,ca_10_loop_1_e39,ca_10_loop_1_e38")
call coverLoop(oEditor, "ca_10_loop_1_e110")
call rename(oEditor, "ca_10_loop_1_e110", "ca_10")
call changeObjectColor(oEditor, "ca_10", 255,127,39)
call drawSegment(oEditor, 42, 45, "ca_11_loop_1_e111")
call drawSegment(oEditor, 44, 45, "ca_11_loop_1_e44")
call drawSegment(oEditor, 43, 44, "ca_11_loop_1_e43")
call drawSegment(oEditor, 42, 43, "ca_11_loop_1_e42")
call uniteEdges(oEditor, "ca_11_loop_1_e111,ca_11_loop_1_e44,ca_11_loop_1_e43,ca_11_loop_1_e42")
call coverLoop(oEditor, "ca_11_loop_1_e111")
call rename(oEditor, "ca_11_loop_1_e111", "ca_11")
call changeObjectColor(oEditor, "ca_11", 255,127,39)
call drawSegment(oEditor, 46, 49, "ca_12_loop_1_e112")
call drawSegment(oEditor, 48, 49, "ca_12_loop_1_e48")
call drawSegment(oEditor, 47, 48, "ca_12_loop_1_e47")
call drawSegment(oEditor, 46, 47, "ca_12_loop_1_e46")
call uniteEdges(oEditor, "ca_12_loop_1_e112,ca_12_loop_1_e48,ca_12_loop_1_e47,ca_12_loop_1_e46")
call coverLoop(oEditor, "ca_12_loop_1_e112")
call rename(oEditor, "ca_12_loop_1_e112", "ca_12")
call changeObjectColor(oEditor, "ca_12", 255,127,39)
call drawSegment(oEditor, 50, 53, "ca_13_loop_1_e113")
call drawSegment(oEditor, 52, 53, "ca_13_loop_1_e52")
call drawSegment(oEditor, 51, 52, "ca_13_loop_1_e51")
call drawSegment(oEditor, 50, 51, "ca_13_loop_1_e50")
call uniteEdges(oEditor, "ca_13_loop_1_e113,ca_13_loop_1_e52,ca_13_loop_1_e51,ca_13_loop_1_e50")
call coverLoop(oEditor, "ca_13_loop_1_e113")
call rename(oEditor, "ca_13_loop_1_e113", "ca_13")
call changeObjectColor(oEditor, "ca_13", 255,127,39)
call drawSegment(oEditor, 54, 57, "ca_14_loop_1_e114")
call drawSegment(oEditor, 56, 57, "ca_14_loop_1_e56")
call drawSegment(oEditor, 55, 56, "ca_14_loop_1_e55")
call drawSegment(oEditor, 54, 55, "ca_14_loop_1_e54")
call uniteEdges(oEditor, "ca_14_loop_1_e114,ca_14_loop_1_e56,ca_14_loop_1_e55,ca_14_loop_1_e54")
call coverLoop(oEditor, "ca_14_loop_1_e114")
call rename(oEditor, "ca_14_loop_1_e114", "ca_14")
call changeObjectColor(oEditor, "ca_14", 255,127,39)
call drawSegment(oEditor, 58, 61, "ca_15_loop_1_e115")
call drawSegment(oEditor, 60, 61, "ca_15_loop_1_e60")
call drawSegment(oEditor, 59, 60, "ca_15_loop_1_e59")
call drawSegment(oEditor, 58, 59, "ca_15_loop_1_e58")
call uniteEdges(oEditor, "ca_15_loop_1_e115,ca_15_loop_1_e60,ca_15_loop_1_e59,ca_15_loop_1_e58")
call coverLoop(oEditor, "ca_15_loop_1_e115")
call rename(oEditor, "ca_15_loop_1_e115", "ca_15")
call changeObjectColor(oEditor, "ca_15", 255,127,39)
call drawSegment(oEditor, 62, 65, "ca_16_loop_1_e116")
call drawSegment(oEditor, 64, 65, "ca_16_loop_1_e64")
call drawSegment(oEditor, 63, 64, "ca_16_loop_1_e63")
call drawSegment(oEditor, 62, 63, "ca_16_loop_1_e62")
call uniteEdges(oEditor, "ca_16_loop_1_e116,ca_16_loop_1_e64,ca_16_loop_1_e63,ca_16_loop_1_e62")
call coverLoop(oEditor, "ca_16_loop_1_e116")
call rename(oEditor, "ca_16_loop_1_e116", "ca_16")
call changeObjectColor(oEditor, "ca_16", 255,127,39)
call drawSegment(oEditor, 66, 69, "ca_17_loop_1_e117")
call drawSegment(oEditor, 68, 69, "ca_17_loop_1_e68")
call drawSegment(oEditor, 67, 68, "ca_17_loop_1_e67")
call drawSegment(oEditor, 66, 67, "ca_17_loop_1_e66")
call uniteEdges(oEditor, "ca_17_loop_1_e117,ca_17_loop_1_e68,ca_17_loop_1_e67,ca_17_loop_1_e66")
call coverLoop(oEditor, "ca_17_loop_1_e117")
call rename(oEditor, "ca_17_loop_1_e117", "ca_17")
call changeObjectColor(oEditor, "ca_17", 255,127,39)
call drawSegment(oEditor, 70, 73, "ca_18_loop_1_e118")
call drawSegment(oEditor, 72, 73, "ca_18_loop_1_e72")
call drawSegment(oEditor, 71, 72, "ca_18_loop_1_e71")
call drawSegment(oEditor, 70, 71, "ca_18_loop_1_e70")
call uniteEdges(oEditor, "ca_18_loop_1_e118,ca_18_loop_1_e72,ca_18_loop_1_e71,ca_18_loop_1_e70")
call coverLoop(oEditor, "ca_18_loop_1_e118")
call rename(oEditor, "ca_18_loop_1_e118", "ca_18")
call changeObjectColor(oEditor, "ca_18", 255,127,39)
call drawSegment(oEditor, 74, 77, "ca_19_loop_1_e119")
call drawSegment(oEditor, 76, 77, "ca_19_loop_1_e76")
call drawSegment(oEditor, 75, 76, "ca_19_loop_1_e75")
call drawSegment(oEditor, 74, 75, "ca_19_loop_1_e74")
call uniteEdges(oEditor, "ca_19_loop_1_e119,ca_19_loop_1_e76,ca_19_loop_1_e75,ca_19_loop_1_e74")
call coverLoop(oEditor, "ca_19_loop_1_e119")
call rename(oEditor, "ca_19_loop_1_e119", "ca_19")
call changeObjectColor(oEditor, "ca_19", 255,127,39)
call drawSegment(oEditor, 78, 81, "ca_20_loop_1_e120")
call drawSegment(oEditor, 80, 81, "ca_20_loop_1_e80")
call drawSegment(oEditor, 79, 80, "ca_20_loop_1_e79")
call drawSegment(oEditor, 78, 79, "ca_20_loop_1_e78")
call uniteEdges(oEditor, "ca_20_loop_1_e120,ca_20_loop_1_e80,ca_20_loop_1_e79,ca_20_loop_1_e78")
call coverLoop(oEditor, "ca_20_loop_1_e120")
call rename(oEditor, "ca_20_loop_1_e120", "ca_20")
call changeObjectColor(oEditor, "ca_20", 255,127,39)
call drawSegment(oEditor, 82, 85, "ca_21_loop_1_e121")
call drawSegment(oEditor, 84, 85, "ca_21_loop_1_e84")
call drawSegment(oEditor, 83, 84, "ca_21_loop_1_e83")
call drawSegment(oEditor, 82, 83, "ca_21_loop_1_e82")
call uniteEdges(oEditor, "ca_21_loop_1_e121,ca_21_loop_1_e84,ca_21_loop_1_e83,ca_21_loop_1_e82")
call coverLoop(oEditor, "ca_21_loop_1_e121")
call rename(oEditor, "ca_21_loop_1_e121", "ca_21")
call changeObjectColor(oEditor, "ca_21", 255,127,39)
call drawSegment(oEditor, 86, 89, "ca_22_loop_1_e122")
call drawSegment(oEditor, 88, 89, "ca_22_loop_1_e88")
call drawSegment(oEditor, 87, 88, "ca_22_loop_1_e87")
call drawSegment(oEditor, 86, 87, "ca_22_loop_1_e86")
call uniteEdges(oEditor, "ca_22_loop_1_e122,ca_22_loop_1_e88,ca_22_loop_1_e87,ca_22_loop_1_e86")
call coverLoop(oEditor, "ca_22_loop_1_e122")
call rename(oEditor, "ca_22_loop_1_e122", "ca_22")
call changeObjectColor(oEditor, "ca_22", 255,127,39)
call drawSegment(oEditor, 90, 93, "ca_23_loop_1_e123")
call drawSegment(oEditor, 92, 93, "ca_23_loop_1_e92")
call drawSegment(oEditor, 91, 92, "ca_23_loop_1_e91")
call drawSegment(oEditor, 90, 91, "ca_23_loop_1_e90")
call uniteEdges(oEditor, "ca_23_loop_1_e123,ca_23_loop_1_e92,ca_23_loop_1_e91,ca_23_loop_1_e90")
call coverLoop(oEditor, "ca_23_loop_1_e123")
call rename(oEditor, "ca_23_loop_1_e123", "ca_23")
call changeObjectColor(oEditor, "ca_23", 255,127,39)
call drawSegment(oEditor, 94, 97, "ca_24_loop_1_e124")
call drawSegment(oEditor, 96, 97, "ca_24_loop_1_e96")
call drawSegment(oEditor, 95, 96, "ca_24_loop_1_e95")
call drawSegment(oEditor, 94, 95, "ca_24_loop_1_e94")
call uniteEdges(oEditor, "ca_24_loop_1_e124,ca_24_loop_1_e96,ca_24_loop_1_e95,ca_24_loop_1_e94")
call coverLoop(oEditor, "ca_24_loop_1_e124")
call rename(oEditor, "ca_24_loop_1_e124", "ca_24")
call changeObjectColor(oEditor, "ca_24", 255,127,39)
call drawSegment(oEditor, 98, 101, "moving_air_loop_1_e125")
call drawSegment(oEditor, 101, 102, "moving_air_loop_1_e126")
call drawSegment(oEditor, 102, 103, "moving_air_loop_1_e127")
call drawSegment(oEditor, 103, 104, "moving_air_loop_1_e128")
call drawSegment(oEditor, 104, 1, "moving_air_loop_1_e129")
call drawSegment(oEditor, 100, 1, "moving_air_loop_1_e100")
call drawSegment(oEditor, 99, 100, "moving_air_loop_1_e99")
call drawSegment(oEditor, 98, 99, "moving_air_loop_1_e98")
call uniteEdges(oEditor, "moving_air_loop_1_e125,moving_air_loop_1_e126,moving_air_loop_1_e127,moving_air_loop_1_e128,moving_air_loop_1_e129,moving_air_loop_1_e100,moving_air_loop_1_e99,moving_air_loop_1_e98")
call coverLoop(oEditor, "moving_air_loop_1_e125")
call rename(oEditor, "moving_air_loop_1_e125", "moving_air")
call changeObjectColor(oEditor, "moving_air", 0,255,255)
call drawSegment(oEditor, 1, 2, "moving_core_loop_1_e1")
call drawSegment(oEditor, 2, 3, "moving_core_loop_1_e2")
call drawSegment(oEditor, 3, 4, "moving_core_loop_1_e3")
call drawSegment(oEditor, 4, 5, "moving_core_loop_1_e4")
call drawSegment(oEditor, 5, 6, "moving_core_loop_1_e5")
call drawSegment(oEditor, 6, 7, "moving_core_loop_1_e6")
call drawSegment(oEditor, 7, 8, "moving_core_loop_1_e7")
call drawSegment(oEditor, 8, 9, "moving_core_loop_1_e8")
call drawSegment(oEditor, 9, 10, "moving_core_loop_1_e9")
call drawSegment(oEditor, 10, 11, "moving_core_loop_1_e10")
call drawSegment(oEditor, 11, 12, "moving_core_loop_1_e11")
call drawSegment(oEditor, 12, 13, "moving_core_loop_1_e12")
call drawSegment(oEditor, 13, 14, "moving_core_loop_1_e13")
call drawSegment(oEditor, 14, 15, "moving_core_loop_1_e14")
call drawSegment(oEditor, 15, 16, "moving_core_loop_1_e15")
call drawSegment(oEditor, 16, 17, "moving_core_loop_1_e16")
call drawSegment(oEditor, 17, 18, "moving_core_loop_1_e17")
call drawSegment(oEditor, 18, 19, "moving_core_loop_1_e18")
call drawSegment(oEditor, 19, 20, "moving_core_loop_1_e19")
call drawSegment(oEditor, 20, 21, "moving_core_loop_1_e20")
call drawSegment(oEditor, 21, 22, "moving_core_loop_1_e21")
call drawSegment(oEditor, 22, 23, "moving_core_loop_1_e22")
call drawSegment(oEditor, 23, 24, "moving_core_loop_1_e23")
call drawSegment(oEditor, 24, 25, "moving_core_loop_1_e24")
call drawSegment(oEditor, 25, 26, "moving_core_loop_1_e25")
call drawSegment(oEditor, 26, 27, "moving_core_loop_1_e26")
call drawSegment(oEditor, 27, 28, "moving_core_loop_1_e27")
call drawSegment(oEditor, 28, 29, "moving_core_loop_1_e28")
call drawSegment(oEditor, 29, 30, "moving_core_loop_1_e29")
call drawSegment(oEditor, 30, 31, "moving_core_loop_1_e30")
call drawSegment(oEditor, 31, 32, "moving_core_loop_1_e31")
call drawSegment(oEditor, 32, 33, "moving_core_loop_1_e32")
call drawSegment(oEditor, 33, 34, "moving_core_loop_1_e33")
call drawSegment(oEditor, 34, 35, "moving_core_loop_1_e34")
call drawSegment(oEditor, 35, 36, "moving_core_loop_1_e35")
call drawSegment(oEditor, 36, 37, "moving_core_loop_1_e36")
call drawSegment(oEditor, 37, 38, "moving_core_loop_1_e37")
call drawSegment(oEditor, 38, 39, "moving_core_loop_1_e38")
call drawSegment(oEditor, 39, 40, "moving_core_loop_1_e39")
call drawSegment(oEditor, 40, 41, "moving_core_loop_1_e40")
call drawSegment(oEditor, 41, 42, "moving_core_loop_1_e41")
call drawSegment(oEditor, 42, 43, "moving_core_loop_1_e42")
call drawSegment(oEditor, 43, 44, "moving_core_loop_1_e43")
call drawSegment(oEditor, 44, 45, "moving_core_loop_1_e44")
call drawSegment(oEditor, 45, 46, "moving_core_loop_1_e45")
call drawSegment(oEditor, 46, 47, "moving_core_loop_1_e46")
call drawSegment(oEditor, 47, 48, "moving_core_loop_1_e47")
call drawSegment(oEditor, 48, 49, "moving_core_loop_1_e48")
call drawSegment(oEditor, 49, 50, "moving_core_loop_1_e49")
call drawSegment(oEditor, 50, 51, "moving_core_loop_1_e50")
call drawSegment(oEditor, 51, 52, "moving_core_loop_1_e51")
call drawSegment(oEditor, 52, 53, "moving_core_loop_1_e52")
call drawSegment(oEditor, 53, 54, "moving_core_loop_1_e53")
call drawSegment(oEditor, 54, 55, "moving_core_loop_1_e54")
call drawSegment(oEditor, 55, 56, "moving_core_loop_1_e55")
call drawSegment(oEditor, 56, 57, "moving_core_loop_1_e56")
call drawSegment(oEditor, 57, 58, "moving_core_loop_1_e57")
call drawSegment(oEditor, 58, 59, "moving_core_loop_1_e58")
call drawSegment(oEditor, 59, 60, "moving_core_loop_1_e59")
call drawSegment(oEditor, 60, 61, "moving_core_loop_1_e60")
call drawSegment(oEditor, 61, 62, "moving_core_loop_1_e61")
call drawSegment(oEditor, 62, 63, "moving_core_loop_1_e62")
call drawSegment(oEditor, 63, 64, "moving_core_loop_1_e63")
call drawSegment(oEditor, 64, 65, "moving_core_loop_1_e64")
call drawSegment(oEditor, 65, 66, "moving_core_loop_1_e65")
call drawSegment(oEditor, 66, 67, "moving_core_loop_1_e66")
call drawSegment(oEditor, 67, 68, "moving_core_loop_1_e67")
call drawSegment(oEditor, 68, 69, "moving_core_loop_1_e68")
call drawSegment(oEditor, 69, 70, "moving_core_loop_1_e69")
call drawSegment(oEditor, 70, 71, "moving_core_loop_1_e70")
call drawSegment(oEditor, 71, 72, "moving_core_loop_1_e71")
call drawSegment(oEditor, 72, 73, "moving_core_loop_1_e72")
call drawSegment(oEditor, 73, 74, "moving_core_loop_1_e73")
call drawSegment(oEditor, 74, 75, "moving_core_loop_1_e74")
call drawSegment(oEditor, 75, 76, "moving_core_loop_1_e75")
call drawSegment(oEditor, 76, 77, "moving_core_loop_1_e76")
call drawSegment(oEditor, 77, 78, "moving_core_loop_1_e77")
call drawSegment(oEditor, 78, 79, "moving_core_loop_1_e78")
call drawSegment(oEditor, 79, 80, "moving_core_loop_1_e79")
call drawSegment(oEditor, 80, 81, "moving_core_loop_1_e80")
call drawSegment(oEditor, 81, 82, "moving_core_loop_1_e81")
call drawSegment(oEditor, 82, 83, "moving_core_loop_1_e82")
call drawSegment(oEditor, 83, 84, "moving_core_loop_1_e83")
call drawSegment(oEditor, 84, 85, "moving_core_loop_1_e84")
call drawSegment(oEditor, 85, 86, "moving_core_loop_1_e85")
call drawSegment(oEditor, 86, 87, "moving_core_loop_1_e86")
call drawSegment(oEditor, 87, 88, "moving_core_loop_1_e87")
call drawSegment(oEditor, 88, 89, "moving_core_loop_1_e88")
call drawSegment(oEditor, 89, 90, "moving_core_loop_1_e89")
call drawSegment(oEditor, 90, 91, "moving_core_loop_1_e90")
call drawSegment(oEditor, 91, 92, "moving_core_loop_1_e91")
call drawSegment(oEditor, 92, 93, "moving_core_loop_1_e92")
call drawSegment(oEditor, 93, 94, "moving_core_loop_1_e93")
call drawSegment(oEditor, 94, 95, "moving_core_loop_1_e94")
call drawSegment(oEditor, 95, 96, "moving_core_loop_1_e95")
call drawSegment(oEditor, 96, 97, "moving_core_loop_1_e96")
call drawSegment(oEditor, 97, 98, "moving_core_loop_1_e97")
call drawSegment(oEditor, 98, 99, "moving_core_loop_1_e98")
call drawSegment(oEditor, 99, 100, "moving_core_loop_1_e99")
call drawSegment(oEditor, 100, 1, "moving_core_loop_1_e100")
call uniteEdges(oEditor, "moving_core_loop_1_e1,moving_core_loop_1_e2,moving_core_loop_1_e3,moving_core_loop_1_e4,moving_core_loop_1_e5,moving_core_loop_1_e6,moving_core_loop_1_e7,moving_core_loop_1_e8,moving_core_loop_1_e9,moving_core_loop_1_e10,moving_core_loop_1_e11,moving_core_loop_1_e12,moving_core_loop_1_e13,moving_core_loop_1_e14,moving_core_loop_1_e15,moving_core_loop_1_e16,moving_core_loop_1_e17,moving_core_loop_1_e18,moving_core_loop_1_e19,moving_core_loop_1_e20,moving_core_loop_1_e21,moving_core_loop_1_e22,moving_core_loop_1_e23,moving_core_loop_1_e24,moving_core_loop_1_e25,moving_core_loop_1_e26,moving_core_loop_1_e27,moving_core_loop_1_e28,moving_core_loop_1_e29,moving_core_loop_1_e30,moving_core_loop_1_e31,moving_core_loop_1_e32,moving_core_loop_1_e33,moving_core_loop_1_e34,moving_core_loop_1_e35,moving_core_loop_1_e36,moving_core_loop_1_e37,moving_core_loop_1_e38,moving_core_loop_1_e39,moving_core_loop_1_e40,moving_core_loop_1_e41,moving_core_loop_1_e42,moving_core_loop_1_e43,moving_core_loop_1_e44,moving_core_loop_1_e45,moving_core_loop_1_e46,moving_core_loop_1_e47,moving_core_loop_1_e48,moving_core_loop_1_e49,moving_core_loop_1_e50,moving_core_loop_1_e51,moving_core_loop_1_e52,moving_core_loop_1_e53,moving_core_loop_1_e54,moving_core_loop_1_e55,moving_core_loop_1_e56,moving_core_loop_1_e57,moving_core_loop_1_e58,moving_core_loop_1_e59,moving_core_loop_1_e60,moving_core_loop_1_e61,moving_core_loop_1_e62,moving_core_loop_1_e63,moving_core_loop_1_e64,moving_core_loop_1_e65,moving_core_loop_1_e66,moving_core_loop_1_e67,moving_core_loop_1_e68,moving_core_loop_1_e69,moving_core_loop_1_e70,moving_core_loop_1_e71,moving_core_loop_1_e72,moving_core_loop_1_e73,moving_core_loop_1_e74,moving_core_loop_1_e75,moving_core_loop_1_e76,moving_core_loop_1_e77,moving_core_loop_1_e78,moving_core_loop_1_e79,moving_core_loop_1_e80,moving_core_loop_1_e81,moving_core_loop_1_e82,moving_core_loop_1_e83,moving_core_loop_1_e84,moving_core_loop_1_e85,moving_core_loop_1_e86,moving_core_loop_1_e87,moving_core_loop_1_e88,moving_core_loop_1_e89,moving_core_loop_1_e90,moving_core_loop_1_e91,moving_core_loop_1_e92,moving_core_loop_1_e93,moving_core_loop_1_e94,moving_core_loop_1_e95,moving_core_loop_1_e96,moving_core_loop_1_e97,moving_core_loop_1_e98,moving_core_loop_1_e99,moving_core_loop_1_e100")
call coverLoop(oEditor, "moving_core_loop_1_e1")
call rename(oEditor, "moving_core_loop_1_e1", "moving_core")
call changeObjectColor(oEditor, "moving_core", 200,200,200)
call drawSegment(oEditor, 105, 106, "rail_air_loop_1_e130")
call drawSegment(oEditor, 106, 108, "rail_air_loop_1_e137")
call drawSegment(oEditor, 107, 108, "rail_air_loop_1_e131")
call drawSegment(oEditor, 105, 107, "rail_air_loop_1_e134")
call uniteEdges(oEditor, "rail_air_loop_1_e130,rail_air_loop_1_e137,rail_air_loop_1_e131,rail_air_loop_1_e134")
call coverLoop(oEditor, "rail_air_loop_1_e130")
call rename(oEditor, "rail_air_loop_1_e130", "rail_air")
call changeObjectColor(oEditor, "rail_air", 0,255,255)
call drawSegment(oEditor, 107, 108, "rail_core_loop_1_e131")
call drawSegment(oEditor, 108, 110, "rail_core_loop_1_e138")
call drawSegment(oEditor, 109, 110, "rail_core_loop_1_e132")
call drawSegment(oEditor, 107, 109, "rail_core_loop_1_e135")
call uniteEdges(oEditor, "rail_core_loop_1_e131,rail_core_loop_1_e138,rail_core_loop_1_e132,rail_core_loop_1_e135")
call coverLoop(oEditor, "rail_core_loop_1_e131")
call rename(oEditor, "rail_core_loop_1_e131", "rail_core")
call changeObjectColor(oEditor, "rail_core", 200,200,200)
call drawSegment(oEditor, 109, 110, "rail_caluminum_loop_1_e132")
call drawSegment(oEditor, 110, 112, "rail_caluminum_loop_1_e139")
call drawSegment(oEditor, 111, 112, "rail_caluminum_loop_1_e133")
call drawSegment(oEditor, 109, 111, "rail_caluminum_loop_1_e136")
call uniteEdges(oEditor, "rail_caluminum_loop_1_e132,rail_caluminum_loop_1_e139,rail_caluminum_loop_1_e133,rail_caluminum_loop_1_e136")
call coverLoop(oEditor, "rail_caluminum_loop_1_e132")
call rename(oEditor, "rail_caluminum_loop_1_e132", "rail_caluminum")
call changeObjectColor(oEditor, "rail_caluminum", 37,177,76)

