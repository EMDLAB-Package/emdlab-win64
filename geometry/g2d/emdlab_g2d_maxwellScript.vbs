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

call defineGlobalVariable(oProject, "x_pts", "[0,42.4933818376,43.4932261162,43.7555069467,54.488544792,57.5,70,69.7336288664,42.3382746689,42.5,43.7933614075,54.7299689092,17.5,23.6704985143,27.1704985143,41.5,35.9400542571,15.1554445662,36.9369703348,35.1869703348,26.5267162969,28.2767162969,27.1704985143,23.6704985143,27.1704985143,27.6704985143,25.9204985143,23.6704985143,36.0619703348] mm")
call makeGBHidden(oProject, "x_pts")
call defineGlobalVariable(oProject, "y_pts", "[0,0.75,0.767647058824,1.82047114958,2.75949028624,0,0,6.10090199234,3.70411906678,0,0,0,0,0,0,0,20.75,8.75,12.665315728,15.6964046412,10.6964046412,7.66531572798,5.74929032419,5.74929032419,6.44929032419,7.31531572798,10.3464046412,6.44929032419,14.1808601846] mm")
call makeGBHidden(oProject, "y_pts")
call defineGlobalVariable(oProject, "e_angles", "[0,0,0,-95,0,5,0,-3.98884552189,0,0,1.01115447811,2.38244734047,0,0,0,30,0,-30,0,0,0,0,0,0,0,0,0,0,0,180] deg")
call makeGBHidden(oProject, "e_angles")
call drawSegment(oEditor, 2, 3, "stator_loop_1_e1")
call drawSegment(oEditor, 3, 4, "stator_loop_1_e2")
call drawSegment(oEditor, 4, 5, "stator_loop_1_e3")
call drawArcCPA(oEditor, 12, 5, 4, "stator_loop_1_e4")
call drawSegment(oEditor, 6, 7, "stator_loop_1_e5")
call drawArcCPA(oEditor, 1, 7, 6, "stator_loop_1_e6")
call drawSegment(oEditor, 8, 9, "stator_loop_1_e7")
call drawArcCPA(oEditor, 1, 9, 8, "stator_loop_1_e8")
call uniteEdges(oEditor, "stator_loop_1_e1,stator_loop_1_e2,stator_loop_1_e3,stator_loop_1_e4,stator_loop_1_e5,stator_loop_1_e6,stator_loop_1_e7,stator_loop_1_e8")
call coverLoop(oEditor, "stator_loop_1_e1")
call rename(oEditor, "stator_loop_1_e1", "stator")
call changeObjectColor(oEditor, "stator", 200,200,200)
call drawSegment(oEditor, 11, 6, "sc_loop_1_e10")
call drawArcCPA(oEditor, 12, 5, 4, "sc_loop_1_e4")
call drawSegment(oEditor, 4, 5, "sc_loop_1_e3")
call drawArcCPA(oEditor, 1, 11, 12, "sc_loop_1_e12")
call uniteEdges(oEditor, "sc_loop_1_e10,sc_loop_1_e4,sc_loop_1_e3,sc_loop_1_e12")
call coverLoop(oEditor, "sc_loop_1_e10")
call rename(oEditor, "sc_loop_1_e10", "sc")
call changeObjectColor(oEditor, "sc", 255,137,39)
call drawSegment(oEditor, 10, 11, "sap_loop_1_e9")
call drawArcCPA(oEditor, 1, 11, 12, "sap_loop_1_e12")
call drawSegment(oEditor, 3, 4, "sap_loop_1_e2")
call drawSegment(oEditor, 2, 3, "sap_loop_1_e1")
call drawArcCPA(oEditor, 1, 10, 11, "sap_loop_1_e11")
call uniteEdges(oEditor, "sap_loop_1_e9,sap_loop_1_e12,sap_loop_1_e2,sap_loop_1_e1,sap_loop_1_e11")
call coverLoop(oEditor, "sap_loop_1_e9")
call rename(oEditor, "sap_loop_1_e9", "sap")
call changeObjectColor(oEditor, "sap", 0,255,255)
call drawSegment(oEditor, 13, 14, "rotor_loop_1_e13")
call drawSegment(oEditor, 24, 14, "rotor_loop_1_e21")
call drawSegment(oEditor, 23, 24, "rotor_loop_1_e20")
call drawSegment(oEditor, 15, 23, "rotor_loop_1_e19")
call drawSegment(oEditor, 15, 16, "rotor_loop_1_e15")
call drawArcCPA(oEditor, 1, 16, 16, "rotor_loop_1_e16")
call drawSegment(oEditor, 17, 18, "rotor_loop_1_e17")
call drawArcCPA(oEditor, 1, 18, 18, "rotor_loop_1_e18")
call uniteEdges(oEditor, "rotor_loop_1_e13,rotor_loop_1_e21,rotor_loop_1_e20,rotor_loop_1_e19,rotor_loop_1_e15,rotor_loop_1_e16,rotor_loop_1_e17,rotor_loop_1_e18")
call coverLoop(oEditor, "rotor_loop_1_e13")
call drawSegment(oEditor, 25, 26, "rotor_loop_2_e22")
call drawSegment(oEditor, 26, 27, "rotor_loop_2_e23")
call drawSegment(oEditor, 27, 28, "rotor_loop_2_e24")
call drawSegment(oEditor, 28, 25, "rotor_loop_2_e25")
call uniteEdges(oEditor, "rotor_loop_2_e22,rotor_loop_2_e23,rotor_loop_2_e24,rotor_loop_2_e25")
call coverLoop(oEditor, "rotor_loop_2_e22")
call drawSegment(oEditor, 22, 19, "rotor_loop_3_e29")
call drawArcCPA(oEditor, 29, 19, 30, "rotor_loop_3_e30")
call drawSegment(oEditor, 20, 21, "rotor_loop_3_e27")
call drawSegment(oEditor, 21, 22, "rotor_loop_3_e28")
call uniteEdges(oEditor, "rotor_loop_3_e29,rotor_loop_3_e30,rotor_loop_3_e27,rotor_loop_3_e28")
call coverLoop(oEditor, "rotor_loop_3_e29")
call subtract(oEditor, "rotor_loop_2_e22,rotor_loop_3_e29", "rotor_loop_1_e13")
call rename(oEditor, "rotor_loop_1_e13", "rotor")
call changeObjectColor(oEditor, "rotor", 200,200,200)
call drawSegment(oEditor, 14, 15, "magnet1_loop_1_e14")
call drawSegment(oEditor, 15, 23, "magnet1_loop_1_e19")
call drawSegment(oEditor, 23, 24, "magnet1_loop_1_e20")
call drawSegment(oEditor, 24, 14, "magnet1_loop_1_e21")
call uniteEdges(oEditor, "magnet1_loop_1_e14,magnet1_loop_1_e19,magnet1_loop_1_e20,magnet1_loop_1_e21")
call coverLoop(oEditor, "magnet1_loop_1_e14")
call rename(oEditor, "magnet1_loop_1_e14", "magnet1")
call changeObjectColor(oEditor, "magnet1", 160,78,146)
call drawSegment(oEditor, 19, 20, "magnet2_loop_1_e26")
call drawSegment(oEditor, 20, 21, "magnet2_loop_1_e27")
call drawSegment(oEditor, 21, 22, "magnet2_loop_1_e28")
call drawSegment(oEditor, 22, 19, "magnet2_loop_1_e29")
call uniteEdges(oEditor, "magnet2_loop_1_e26,magnet2_loop_1_e27,magnet2_loop_1_e28,magnet2_loop_1_e29")
call coverLoop(oEditor, "magnet2_loop_1_e26")
call rename(oEditor, "magnet2_loop_1_e26", "magnet2")
call changeObjectColor(oEditor, "magnet2", 160,78,146)
call drawSegment(oEditor, 25, 26, "rap1_loop_1_e22")
call drawSegment(oEditor, 26, 27, "rap1_loop_1_e23")
call drawSegment(oEditor, 27, 28, "rap1_loop_1_e24")
call drawSegment(oEditor, 28, 25, "rap1_loop_1_e25")
call uniteEdges(oEditor, "rap1_loop_1_e22,rap1_loop_1_e23,rap1_loop_1_e24,rap1_loop_1_e25")
call coverLoop(oEditor, "rap1_loop_1_e22")
call rename(oEditor, "rap1_loop_1_e22", "rap1")
call changeObjectColor(oEditor, "rap1", 0,255,255)
call drawArcCPA(oEditor, 29, 19, 30, "rap2_loop_1_e30")
call drawSegment(oEditor, 19, 20, "rap2_loop_1_e26")
call uniteEdges(oEditor, "rap2_loop_1_e30,rap2_loop_1_e26")
call coverLoop(oEditor, "rap2_loop_1_e30")
call rename(oEditor, "rap2_loop_1_e30", "rap2")
call changeObjectColor(oEditor, "rap2", 0,255,255)

