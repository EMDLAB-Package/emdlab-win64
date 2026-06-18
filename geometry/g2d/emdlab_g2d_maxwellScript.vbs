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

call defineGlobalVariable(oProject, "x_pts", "[47.9948955619,48.9948955619,49.484518603,62.784518603,62.784518603,70,0,69.8501246267,47.8972283155,49.984518603,52.684518603,52.684518603,49.984518603,53.184518603,55.884518603,55.884518603,53.184518603,56.384518603,59.084518603,59.084518603,56.384518603,59.584518603,62.284518603,62.284518603,59.584518603,48,32,47.2,44.7962562113,42.683234402,29.5641450404,39.7949731348,43.7949731348,43.7949731348,40.7949731348,39.7949731348,43.7949731348,42.7463888848,42.0662123714,41.4142363709,40.7949731348,43.6071139345] mm")
call makeGBHidden(oProject, "x_pts")
call defineGlobalVariable(oProject, "y_pts", "[0.7,0.7,1.75,1.75,-1.04360964315e-14,0,0,4.57821904611,3.13935020305,0,0,1.55,1.55,0,0,1.55,1.55,0,0,1.55,1.55,0,0,1.55,1.55,0,0,0,14.8706230349,17.6799745753,12.2458698357,0.5,0.5,14.2298493675,14.2298493675,14.2298493675,14.5382358124,15.932372761,15.4207857324,15.9678585539,15.2298493675,18.0626580076] mm")
call makeGBHidden(oProject, "y_pts")
call defineGlobalVariable(oProject, "e_angles", "[0,0,0,0,0,3.75,0,-2.91440692888,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.835593071117,0,0,18.3641687031,0,0,-22.5,0,0,0,0,0,0,0,103.051760884,0,0,4.13583129689,0] deg")
call makeGBHidden(oProject, "e_angles")
call drawSegment(oEditor, 10, 11, "sc1_loop_1_e9")
call drawSegment(oEditor, 11, 12, "sc1_loop_1_e10")
call drawSegment(oEditor, 12, 13, "sc1_loop_1_e11")
call drawSegment(oEditor, 13, 10, "sc1_loop_1_e12")
call uniteEdges(oEditor, "sc1_loop_1_e9,sc1_loop_1_e10,sc1_loop_1_e11,sc1_loop_1_e12")
call coverLoop(oEditor, "sc1_loop_1_e9")
call rename(oEditor, "sc1_loop_1_e9", "sc1")
call changeObjectColor(oEditor, "sc1", 255,137,39)
call drawSegment(oEditor, 14, 15, "sc2_loop_1_e14")
call drawSegment(oEditor, 15, 16, "sc2_loop_1_e15")
call drawSegment(oEditor, 16, 17, "sc2_loop_1_e16")
call drawSegment(oEditor, 17, 14, "sc2_loop_1_e17")
call uniteEdges(oEditor, "sc2_loop_1_e14,sc2_loop_1_e15,sc2_loop_1_e16,sc2_loop_1_e17")
call coverLoop(oEditor, "sc2_loop_1_e14")
call rename(oEditor, "sc2_loop_1_e14", "sc2")
call changeObjectColor(oEditor, "sc2", 255,137,39)
call drawSegment(oEditor, 18, 19, "sc3_loop_1_e19")
call drawSegment(oEditor, 19, 20, "sc3_loop_1_e20")
call drawSegment(oEditor, 20, 21, "sc3_loop_1_e21")
call drawSegment(oEditor, 21, 18, "sc3_loop_1_e22")
call uniteEdges(oEditor, "sc3_loop_1_e19,sc3_loop_1_e20,sc3_loop_1_e21,sc3_loop_1_e22")
call coverLoop(oEditor, "sc3_loop_1_e19")
call rename(oEditor, "sc3_loop_1_e19", "sc3")
call changeObjectColor(oEditor, "sc3", 255,137,39)
call drawSegment(oEditor, 22, 23, "sc4_loop_1_e24")
call drawSegment(oEditor, 23, 24, "sc4_loop_1_e25")
call drawSegment(oEditor, 24, 25, "sc4_loop_1_e26")
call drawSegment(oEditor, 25, 22, "sc4_loop_1_e27")
call uniteEdges(oEditor, "sc4_loop_1_e24,sc4_loop_1_e25,sc4_loop_1_e26,sc4_loop_1_e27")
call coverLoop(oEditor, "sc4_loop_1_e24")
call rename(oEditor, "sc4_loop_1_e24", "sc4")
call changeObjectColor(oEditor, "sc4", 255,137,39)
call drawSegment(oEditor, 1, 2, "stator_loop_1_e1")
call drawSegment(oEditor, 2, 3, "stator_loop_1_e2")
call drawSegment(oEditor, 3, 4, "stator_loop_1_e3")
call drawSegment(oEditor, 4, 5, "stator_loop_1_e4")
call drawSegment(oEditor, 5, 6, "stator_loop_1_e5")
call drawArcCPA(oEditor, 7, 6, 6, "stator_loop_1_e6")
call drawSegment(oEditor, 8, 9, "stator_loop_1_e7")
call drawArcCPA(oEditor, 7, 9, 8, "stator_loop_1_e8")
call uniteEdges(oEditor, "stator_loop_1_e1,stator_loop_1_e2,stator_loop_1_e3,stator_loop_1_e4,stator_loop_1_e5,stator_loop_1_e6,stator_loop_1_e7,stator_loop_1_e8")
call coverLoop(oEditor, "stator_loop_1_e1")
call rename(oEditor, "stator_loop_1_e1", "stator")
call changeObjectColor(oEditor, "stator", 200,200,200)
call drawSegment(oEditor, 13, 10, "sap_loop_1_e12")
call drawSegment(oEditor, 12, 13, "sap_loop_1_e11")
call drawSegment(oEditor, 11, 12, "sap_loop_1_e10")
call drawSegment(oEditor, 11, 14, "sap_loop_1_e13")
call drawSegment(oEditor, 17, 14, "sap_loop_1_e17")
call drawSegment(oEditor, 16, 17, "sap_loop_1_e16")
call drawSegment(oEditor, 15, 16, "sap_loop_1_e15")
call drawSegment(oEditor, 15, 18, "sap_loop_1_e18")
call drawSegment(oEditor, 21, 18, "sap_loop_1_e22")
call drawSegment(oEditor, 20, 21, "sap_loop_1_e21")
call drawSegment(oEditor, 19, 20, "sap_loop_1_e20")
call drawSegment(oEditor, 19, 22, "sap_loop_1_e23")
call drawSegment(oEditor, 25, 22, "sap_loop_1_e27")
call drawSegment(oEditor, 24, 25, "sap_loop_1_e26")
call drawSegment(oEditor, 23, 24, "sap_loop_1_e25")
call drawSegment(oEditor, 23, 5, "sap_loop_1_e28")
call drawSegment(oEditor, 4, 5, "sap_loop_1_e4")
call drawSegment(oEditor, 3, 4, "sap_loop_1_e3")
call drawSegment(oEditor, 2, 3, "sap_loop_1_e2")
call drawSegment(oEditor, 1, 2, "sap_loop_1_e1")
call drawArcCPA(oEditor, 7, 1, 29, "sap_loop_1_e29")
call drawSegment(oEditor, 26, 10, "sap_loop_1_e30")
call uniteEdges(oEditor, "sap_loop_1_e12,sap_loop_1_e11,sap_loop_1_e10,sap_loop_1_e13,sap_loop_1_e17,sap_loop_1_e16,sap_loop_1_e15,sap_loop_1_e18,sap_loop_1_e22,sap_loop_1_e21,sap_loop_1_e20,sap_loop_1_e23,sap_loop_1_e27,sap_loop_1_e26,sap_loop_1_e25,sap_loop_1_e28,sap_loop_1_e4,sap_loop_1_e3,sap_loop_1_e2,sap_loop_1_e1,sap_loop_1_e29,sap_loop_1_e30")
call coverLoop(oEditor, "sap_loop_1_e12")
call rename(oEditor, "sap_loop_1_e12", "sap")
call changeObjectColor(oEditor, "sap", 0,255,255)
call drawSegment(oEditor, 27, 28, "rotor_loop_1_e31")
call drawArcCPA(oEditor, 7, 28, 32, "rotor_loop_1_e32")
call drawSegment(oEditor, 29, 30, "rotor_loop_1_e33")
call drawSegment(oEditor, 30, 31, "rotor_loop_1_e34")
call drawArcCPA(oEditor, 7, 31, 35, "rotor_loop_1_e35")
call uniteEdges(oEditor, "rotor_loop_1_e31,rotor_loop_1_e32,rotor_loop_1_e33,rotor_loop_1_e34,rotor_loop_1_e35")
call coverLoop(oEditor, "rotor_loop_1_e31")
call drawSegment(oEditor, 32, 33, "rotor_loop_2_e36")
call drawSegment(oEditor, 33, 34, "rotor_loop_2_e37")
call drawSegment(oEditor, 34, 37, "rotor_loop_2_e41")
call drawSegment(oEditor, 37, 38, "rotor_loop_2_e42")
call drawArcCPA(oEditor, 39, 38, 43, "rotor_loop_2_e43")
call drawSegment(oEditor, 40, 41, "rotor_loop_2_e44")
call drawSegment(oEditor, 41, 35, "rotor_loop_2_e45")
call drawSegment(oEditor, 35, 36, "rotor_loop_2_e39")
call drawSegment(oEditor, 36, 32, "rotor_loop_2_e40")
call uniteEdges(oEditor, "rotor_loop_2_e36,rotor_loop_2_e37,rotor_loop_2_e41,rotor_loop_2_e42,rotor_loop_2_e43,rotor_loop_2_e44,rotor_loop_2_e45,rotor_loop_2_e39,rotor_loop_2_e40")
call coverLoop(oEditor, "rotor_loop_2_e36")
call subtract(oEditor, "rotor_loop_2_e36", "rotor_loop_1_e31")
call rename(oEditor, "rotor_loop_1_e31", "rotor")
call changeObjectColor(oEditor, "rotor", 200,200,200)
call drawSegment(oEditor, 32, 33, "magnet_loop_1_e36")
call drawSegment(oEditor, 33, 34, "magnet_loop_1_e37")
call drawSegment(oEditor, 34, 35, "magnet_loop_1_e38")
call drawSegment(oEditor, 35, 36, "magnet_loop_1_e39")
call drawSegment(oEditor, 36, 32, "magnet_loop_1_e40")
call uniteEdges(oEditor, "magnet_loop_1_e36,magnet_loop_1_e37,magnet_loop_1_e38,magnet_loop_1_e39,magnet_loop_1_e40")
call coverLoop(oEditor, "magnet_loop_1_e36")
call rename(oEditor, "magnet_loop_1_e36", "magnet")
call changeObjectColor(oEditor, "magnet", 255,137,39)
call drawSegment(oEditor, 34, 37, "rap1_loop_1_e41")
call drawSegment(oEditor, 37, 38, "rap1_loop_1_e42")
call drawArcCPA(oEditor, 39, 38, 43, "rap1_loop_1_e43")
call drawSegment(oEditor, 40, 41, "rap1_loop_1_e44")
call drawSegment(oEditor, 41, 35, "rap1_loop_1_e45")
call drawSegment(oEditor, 34, 35, "rap1_loop_1_e38")
call uniteEdges(oEditor, "rap1_loop_1_e41,rap1_loop_1_e42,rap1_loop_1_e43,rap1_loop_1_e44,rap1_loop_1_e45,rap1_loop_1_e38")
call coverLoop(oEditor, "rap1_loop_1_e41")
call rename(oEditor, "rap1_loop_1_e41", "rap1")
call changeObjectColor(oEditor, "rap1", 0,255,255)
call drawArcCPA(oEditor, 7, 29, 46, "rap2_loop_1_e46")
call drawSegment(oEditor, 42, 30, "rap2_loop_1_e47")
call drawSegment(oEditor, 29, 30, "rap2_loop_1_e33")
call uniteEdges(oEditor, "rap2_loop_1_e46,rap2_loop_1_e47,rap2_loop_1_e33")
call coverLoop(oEditor, "rap2_loop_1_e46")
call rename(oEditor, "rap2_loop_1_e46", "rap2")
call changeObjectColor(oEditor, "rap2", 0,255,255)

