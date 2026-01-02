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

call defineGlobalVariable(oProject, "x_pts", "[0,80.9938269253,81.9937507145,82.0355140669,110.344043268,114,135,134.710954637,80.8265727823,81,82.0505119271,110.568475695,64.8749122227,71.7920677438,77.8082533122,72.3936863007,66.2856842656,71.7002512771,74.7269399062,55,80.25,74.141332484,50.8133742881] mm")
call makeGBHidden(oProject, "x_pts")
call defineGlobalVariable(oProject, "y_pts", "[0,1,1.01234567901,1.56873811671,3.42417714774,0,0,8.82942244607,5.29765346764,0,0,0,0.5,0.5,17.0293340004,2.15293340004,4.37606433166,19.2524649321,27.5682235962,0,0,30.7103454473,21.0475887801] mm")
call makeGBHidden(oProject, "y_pts")
call defineGlobalVariable(oProject, "e_angles", "[0,0,0,-93.75,0,3.75,0,-3.04262672745,0,0,0.707373272549,1.09551479932,0,0,0,0,0,7.9047437659,0,0,0,0,22.5,0,-22.5] deg")
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
call drawSegment(oEditor, 11, 6, "sca_loop_1_e10")
call drawArcCPA(oEditor, 12, 5, 4, "sca_loop_1_e4")
call drawSegment(oEditor, 4, 5, "sca_loop_1_e3")
call drawArcCPA(oEditor, 1, 11, 12, "sca_loop_1_e12")
call uniteEdges(oEditor, "sca_loop_1_e10,sca_loop_1_e4,sca_loop_1_e3,sca_loop_1_e12")
call coverLoop(oEditor, "sca_loop_1_e10")
call rename(oEditor, "sca_loop_1_e10", "sca")
call changeObjectColor(oEditor, "sca", 255,137,39)
call drawSegment(oEditor, 10, 11, "sap_loop_1_e9")
call drawArcCPA(oEditor, 1, 11, 12, "sap_loop_1_e12")
call drawSegment(oEditor, 3, 4, "sap_loop_1_e2")
call drawSegment(oEditor, 2, 3, "sap_loop_1_e1")
call drawArcCPA(oEditor, 1, 10, 11, "sap_loop_1_e11")
call uniteEdges(oEditor, "sap_loop_1_e9,sap_loop_1_e12,sap_loop_1_e2,sap_loop_1_e1,sap_loop_1_e11")
call coverLoop(oEditor, "sap_loop_1_e9")
call rename(oEditor, "sap_loop_1_e9", "sap")
call changeObjectColor(oEditor, "sap", 0,255,255)
call drawSegment(oEditor, 20, 21, "rotor_loop_1_e22")
call drawArcCPA(oEditor, 1, 21, 23, "rotor_loop_1_e23")
call drawSegment(oEditor, 22, 23, "rotor_loop_1_e24")
call drawArcCPA(oEditor, 1, 23, 25, "rotor_loop_1_e25")
call uniteEdges(oEditor, "rotor_loop_1_e22,rotor_loop_1_e23,rotor_loop_1_e24,rotor_loop_1_e25")
call coverLoop(oEditor, "rotor_loop_1_e22")
call drawSegment(oEditor, 13, 14, "rotor_loop_2_e13")
call drawSegment(oEditor, 14, 16, "rotor_loop_2_e14")
call drawSegment(oEditor, 16, 15, "rotor_loop_2_e15")
call drawArcCPA(oEditor, 1, 15, 18, "rotor_loop_2_e18")
call drawSegment(oEditor, 19, 18, "rotor_loop_2_e19")
call drawSegment(oEditor, 18, 17, "rotor_loop_2_e20")
call drawSegment(oEditor, 17, 13, "rotor_loop_2_e21")
call uniteEdges(oEditor, "rotor_loop_2_e13,rotor_loop_2_e14,rotor_loop_2_e15,rotor_loop_2_e18,rotor_loop_2_e19,rotor_loop_2_e20,rotor_loop_2_e21")
call coverLoop(oEditor, "rotor_loop_2_e13")
call subtract(oEditor, "rotor_loop_2_e13", "rotor_loop_1_e22")
call rename(oEditor, "rotor_loop_1_e22", "rotor")
call changeObjectColor(oEditor, "rotor", 200,200,200)
call drawSegment(oEditor, 16, 15, "magnet_loop_1_e15")
call drawSegment(oEditor, 15, 18, "magnet_loop_1_e17")
call drawSegment(oEditor, 18, 17, "magnet_loop_1_e20")
call drawSegment(oEditor, 16, 17, "magnet_loop_1_e16")
call uniteEdges(oEditor, "magnet_loop_1_e15,magnet_loop_1_e17,magnet_loop_1_e20,magnet_loop_1_e16")
call coverLoop(oEditor, "magnet_loop_1_e15")
call rename(oEditor, "magnet_loop_1_e15", "magnet")
call changeObjectColor(oEditor, "magnet", 160,78,146)
call drawArcCPA(oEditor, 1, 15, 18, "rap1_loop_1_e18")
call drawSegment(oEditor, 19, 18, "rap1_loop_1_e19")
call drawSegment(oEditor, 15, 18, "rap1_loop_1_e17")
call uniteEdges(oEditor, "rap1_loop_1_e18,rap1_loop_1_e19,rap1_loop_1_e17")
call coverLoop(oEditor, "rap1_loop_1_e18")
call rename(oEditor, "rap1_loop_1_e18", "rap1")
call changeObjectColor(oEditor, "rap1", 0,255,255)
call drawSegment(oEditor, 13, 14, "rap2_loop_1_e13")
call drawSegment(oEditor, 14, 16, "rap2_loop_1_e14")
call drawSegment(oEditor, 16, 17, "rap2_loop_1_e16")
call drawSegment(oEditor, 17, 13, "rap2_loop_1_e21")
call uniteEdges(oEditor, "rap2_loop_1_e13,rap2_loop_1_e14,rap2_loop_1_e16,rap2_loop_1_e21")
call coverLoop(oEditor, "rap2_loop_1_e13")
call rename(oEditor, "rap2_loop_1_e13", "rap2")
call changeObjectColor(oEditor, "rap2", 0,255,255)

