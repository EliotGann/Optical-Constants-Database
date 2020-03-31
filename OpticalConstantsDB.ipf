#pragma rtGlobals=1	// Use modern global access method.
Menu "Macros"
	submenu "Optical Constants database"
		"Upload NEXAFS to Optical Constant DB", /Q, Execute/P/Q "OpenLoadOCPanel()"
			help={"Load calculated Delta and Betas from a graph into a database of Optical Constants"}
		"Load Optical Constants DB", /Q, Execute/P/Q "LoadOC()"
			help={"Loads opticalconstants.txt from User Procedures"}
	end
End
function loadOC()
	string foldersave = getdatafolder(1)
	setdatafolder root:
	newdatafolder /s/o opticalconstants
	//newpath /m="Select location of optical constant database: opticalconstants.txt " /O /Q opticalconstants
	NewPath/C/Z/Q/O opticalconstants SpecialDirPath("Igor Pro User files",0,0,0)+"User Procedures:Optical Constants"
	GetFileFolderInfo /P=opticalconstants /Q /Z "opticalconstants.txt"
	if(v_flag)
		beep
		doalert 0, "No \"Optical Constants/opticalconstants.txt\" file found!"
	else
		LoadWave/J/D/W/O/Q/G/K=0/A /p=opticalconstants "opticalconstants.txt"
	endif
	setdatafolder foldersave
end	
function getdelta(materialname,wavelength,[pathname,aligned])
	string materialname,pathname
	variable wavelength,aligned // aligned is
	if(paramisdefault(pathname))
		pathname = "root:OpticalConstants"
	endif
	
	string foldersave = getdatafolder(1)
	setdatafolder pathname
	aligned = paramisdefault(aligned) ? -1 : aligned
	if(aligned == 0)
		wave deltaw = $("delta_"+materialname +"perp")
		if(!waveexists(deltaw))
			wave deltaw = $("delta_"+materialname)
		endif
		wave enw = $("energy_"+materialname)
	elseif(aligned == 1)
		wave deltaw = $("delta_"+materialname +"para")
		if(!waveexists(deltaw))
			wave deltaw = $("delta_"+materialname)
		endif
		wave enw = $("energy_"+materialname)
	else
		wave deltaw = $("delta_"+materialname )
		wave enw = $("energy_"+materialname)
	endif
	setdatafolder foldersave
	return interp(1239.84197/wavelength,enw,deltaw)
end
function getbeta(materialname,wavelength,[pathname,aligned])
	string materialname,pathname
	variable wavelength,aligned // aligned is
	if(paramisdefault(pathname))
		pathname = "root:OpticalConstants"
	endif
	
	string foldersave = getdatafolder(1)
	setdatafolder pathname
	aligned = paramisdefault(aligned) ? -1 : aligned
	if(aligned == 0)
		wave/z betaw = $("beta_"+materialname +"perp")
		if(!waveexists(betaw))
			wave betaw = $("beta_"+materialname )
		endif
		wave enw = $("energy_"+materialname)
	elseif(aligned == 1)
		wave/z betaw = $("beta_"+materialname +"para")
		if(!waveexists(betaw))
			wave betaw = $("beta_"+materialname )
		endif
		wave enw = $("energy_"+materialname)
	else
		wave betaw = $("beta_"+materialname )
		wave enw = $("energy_"+materialname )
	endif
	setdatafolder foldersave
	return interp(1239.84197/wavelength,enw,betaw)
end

function contrast(mat1,mat2,en)
	string mat1, mat2
	variable en
	variable d1,b1,d2,b2
	string foldersave = getdatafolder(1)
	setdatafolder root:opticalconstants
	wave beta1 = $("beta_"+ mat1 )
	wave delta1 = $("delta_"+ mat1 )
	wave en1 = $("energy_"+ mat1 )
	wave beta2 = $("beta_"+ mat2 )
	wave delta2 = $("delta_"+ mat2 )
	wave en2 = $("energy_"+ mat2 )
	d1 = interp(en,en1,delta1)
	b1 = interp(en,en1,beta1)
	d2 = interp(en,en2,delta2)
	b2 = interp(en,en2,beta2)
	setdatafolder foldersave
	return ((d1-d2)^2 + (b1-b2)^2) * en^-4
end

function OCsexist(material)
	string material
	string foldersave = getdatafolder(1)
	setdatafolder root:opticalconstants
	variable returnv
	wave/z betapa = $("beta_"+material +"perp")
	wave/z betape = $("beta_"+material +"para")
	wave/z enw = $("energy_"+material )
	wave/z betaw = $("beta_"+material )
	wave/z deltaw = $("delta_"+material )
	if(waveexists(enw) && waveexists(deltaw) && waveexists(betaw))
		if(waveexists(betapa)&&waveexists(betape))
			//aligned optical constants probably exists (at least the betas do)
			returnv = 2
		else
			returnv =  1//basic optical constants exists
		endif
	else
		returnv= 0
	endif
	setdatafolder foldersave
	return returnv
end
function /s listOCs()
	string foldersave = getdatafolder(1)
	if(!datafolderexists("root:opticalconstants"))
		loadOC()
	endif
	setdatafolder root:opticalconstants
	string listofwavesd = wavelist("delta_*",";","")
	string listofwavesb = wavelist("beta_*",";","")
	string listofwavese = wavelist("energy_*",";","")
	setdatafolder foldersave
	variable i
	string materiald, material,listofOCs=""
	for(i=0;i<itemsinlist(listofwavesd);i+=1)
		materiald = stringfromlist(i,listofwavesd)
		material = replacestring("delta_",materiald,"")
		if(findlistitem("beta_"+material, listofwavesb)>=0 &&findlistitem("energy_"+material,listofwavese)>=0 && !stringmatch(material,"*para")&& !stringmatch(material,"*perp") )
			listofOCs = listofOCs + material+";"
		endif
	endfor
	return listofOCs
end

function /S LoadOC_ListTraces()
	string basiclist =TraceNameList("",";",1)
	basiclist = "none;"+basiclist
	return basiclist
end

Function MakeGraph_but(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			
			controlinfo Unaligned_deltaPop
			string unaligneddeltatrace = S_Value
			controlinfo Aligned_deltapop
			string aligneddeltatrace = S_Value
			controlinfo Unaligned_betaPop
			string unalignedbetatrace = S_Value
			controlinfo Aligned_betapop
			string alignedbetatrace = S_Value
			wave o_dunalignedw = tracenametowaveref("",unaligneddeltatrace)
			wave/z o_dunalignedwx = xwavereffromtrace("",unaligneddeltatrace)
			wave o_dalignedw = tracenametowaveref("",aligneddeltatrace)
			wave/z o_dalignedwx = xwavereffromtrace("",aligneddeltatrace)
			wave o_bunalignedw = tracenametowaveref("",unalignedbetatrace)
			wave/z o_bunalignedwx = xwavereffromtrace("",unalignedbetatrace)
			wave o_balignedw = tracenametowaveref("",alignedbetatrace)
			wave/z o_balignedwx = xwavereffromtrace("",alignedbetatrace)
			setdatafolder root:
			newdatafolder /o/s OpticalConstants
			newdatafolder /o/s Calcs
			string /g name=""
			duplicate /o o_dunalignedw, wdcenter,wbcenter,wdaligned,wbaligned,wdalignedpara,wbalignedpara,wddif,wbdif
			if(waveexists(o_dunalignedwx))
				duplicate /o o_dunalignedwx, xaxis
			else
				duplicate /o o_dunalignedw, xaxis
				xaxis=x
			endif
			if(waveexists(o_bunalignedwx))
				wbcenter = interp(xaxis,o_bunalignedwx,o_bunalignedw)	
			else
				wbcenter = o_bunalignedw(x)
			endif
			if(waveexists(o_dalignedwx))
				wdaligned = interp(xaxis,o_dalignedwx,o_dalignedw)
				wdalignedpara = interp(xaxis,o_dalignedwx,o_dalignedw)	
			else
				wdaligned = o_dalignedw(x)
				wdalignedpara = o_dalignedw(x)
			endif
			if(waveexists(o_balignedwx))
				wbaligned = interp(xaxis,o_balignedwx,o_balignedw)
				wbalignedpara = interp(xaxis,o_balignedwx,o_balignedw)
			else
				wbaligned = o_balignedw(x)
				wbalignedpara = o_balignedw(x)
			endif
			wddif = wdcenter-wdaligned
			wbdif = wbcenter-wbaligned
			variable /g calc=1
			variable /g centercalc=1
			dowindow /k OCCalcGraph
			display /k=1 /n=OCCalcGraph wdaligned, wdalignedpara, wdcenter, wbaligned, wbalignedpara, wbcenter vs xaxis
			ModifyGraph /w=OCCalcGraph lsize=1.5, rgb(wdcenter)=(0,0,0), rgb(wbcenter)=(0,0,0), rgb(wdalignedpara)=(0,0,40000), rgb(wbalignedpara)=(0,0,40000)
			ModifyGraph /w=OCCalcGraph lstyle(wdaligned)=3, lstyle(wdcenter)=3,lstyle(wdalignedpara)=3
			SetAxis/A=2 left
			SetAxis bottom 100,600
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SaveOC_button(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			string foldersave = getdatafolder(1)
			svar name = root:OpticalConstants:Calcs:name
			if(strlen(name)<1)
				doprompt "Please enter a name for this material"
				setdatafolder foldersave
				return 0
			endif
			dowindow /F OCCalcGraph
			if(v_flag)
				// there is a graph, use those waves graphed
				setdatafolder root:OpticalConstants:Calcs
				sparcewaves(1e-6)
				wave wdcenter0,wbcenter0,wdaligned0,wbaligned0,wdalignedpara0,wbalignedpara0, xwaveout
				setdatafolder root:opticalconstants
				duplicate /o xwaveout, $("energy_"+name)
				duplicate /o wbcenter0, $("beta_"+name)
				duplicate /o wdcenter0, $("delta_"+name)
				duplicate /o wbalignedpara0, $("beta_"+name+"para")
				duplicate /o wdalignedpara0, $("delta_"+name+"para")
				duplicate /o wbaligned0, $("beta_"+name+"perp")
				duplicate /o wdaligned0, $("delta_"+name+"perp")
			else
				// use the selected waves (unaligned)
				controlinfo Unaligned_deltaPop
				string unaligneddeltatrace = S_Value
				controlinfo Unaligned_betaPop
				string unalignedbetatrace = S_Value
				
				wave/z deltaw = tracenametowaveref("",unaligneddeltatrace)
				wave/z deltawx = xwavereffromtrace("",unaligneddeltatrace)
				wave/z betaw = tracenametowaveref("",unalignedbetatrace)
				wave/z betawx = xwavereffromtrace("",unalignedbetatrace)
				if(!waveexists(deltaw) || ! waveexists(betaw))
					beep
					doalert 0, "Choose valid waves"
					setdatafolder foldersave
					return 0
				endif
				setdatafolder root:opticalconstants
				newdatafolder /o/s Calcs
				duplicate/o deltaw, deltawave, betawave
				wave deltawave
				if(waveexists(deltawx))
					duplicate /o deltawx, deltaxwave
				else
					duplicate/o deltaw, deltaxwave
					deltaxwave = x
				endif
				if(waveexists(betawx))
					betawave = interp(deltaxwave,betawx,betaw)
				else
					duplicate /free betaw, betaoldx
					betaoldx=x
					betawave = interp(deltaxwave,betaoldx,betaw) 
				endif
				display /k=1 /n=OCCalcGraph betawave, deltawave vs deltaxwave
				SetAxis/A=2 left
				SetAxis bottom 250,350
				sparcewaves(1e-6)
				wave xwaveout, betawave0, deltawave0
				setdatafolder root:opticalconstants
				
				duplicate /o xwaveout, $("energy_"+name)
				duplicate /o betawave0, $("beta_"+name)
				duplicate /o deltawave0, $("delta_"+name)
			endif
			dowindow /k OCCalcGraph
			doupdate
			killdatafolder /z Calcs
			newdatafolder /o/s Calcs
			string /g name=""
			setdatafolder ::
			// update the opticalconstants DB file
			NewPath/C/Z/Q/O opticalconstants SpecialDirPath("Igor Pro User files",0,0,0)+"User Procedures:Optical Constants"
			
			Save/O/G/W/U={0,0,1,0}/B /p=opticalconstants wavelist("*",";","") as "OpticalConstants.txt"
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function AdjustOC_slider(sa) : SliderControl
	STRUCT WMSliderAction &sa

	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			nvar/z calc = root:OpticalConstants:Calcs:calc
			nvar/z centercalc = root:OpticalConstants:Calcs:CenterCalc
			wave/z wdc = root:OpticalConstants:Calcs:wdcenter
			wave/z wda = root:OpticalConstants:Calcs:wdaligned
			wave/z wdap = root:OpticalConstants:Calcs:wdalignedpara
			wave/z wdd = root:OpticalConstants:Calcs:wddif
			wave/z wbc = root:OpticalConstants:Calcs:wbcenter
			wave/z wba = root:OpticalConstants:Calcs:wbaligned
			wave/z wbap = root:OpticalConstants:Calcs:wbalignedpara
			wave/z wbd = root:OpticalConstants:Calcs:wbdif
			if(!nvar_exists(Calc) || !nvar_exists(centercalc) || !waveexists(wdc) || !waveexists(wda) || !waveexists(wdd) || !waveexists(wbc) || !waveexists(wba) || !waveexists(wbd) )
				slider $sa.ctrlName value=0
				break
			endif
			calc = (4*sa.curval^3)*centercalc//sign(sa.curval)*
			if( sa.eventCode & 1 ) // value set
				//calculate an appropriate difference spectrum based on the slider position
				wba = wbc+(calc*wbd)
				wda = wdc+(calc*wdd)
				wbap = wbc-(2*calc*wbd)
				wdap = wdc-(2*calc*wdd)
			endif
			if( sa.eventCode & 4 ) // mouse up
				//set the slider value back to 0
				centercalc = calc
				slider $sa.ctrlName value=0
			endif
			break
	endswitch

	return 0
End

function OpenLoadOCPanel()
	string foldersave = getdatafolder(1)
	setdatafolder root:
	newdatafolder /o/s OpticalConstants
	newdatafolder /o/s Calcs
	variable /g calc=1
	variable /g centercalc=1
	string /g name=""
	dowindow /k LoadOCPanel
	NewPanel /k=1 /W=(820,201,1088,453) /n=LoadOCPanel as "Load OC"
	SetDrawLayer UserBack
	DrawText 6,211,"Without graphing only unaligned OC are saved"
	PopupMenu Unaligned_betaPop,pos={28,9},size={235,21},bodyWidth=159,title="unaligned beta:"
	PopupMenu Unaligned_betaPop,mode=1,popvalue="none",value= #"LoadOC_ListTraces()"
	PopupMenu Unaligned_deltaPop,pos={26,34},size={237,21},bodyWidth=159,title="unaligned delta:"
	PopupMenu Unaligned_deltaPop,mode=1,popvalue="none",value= #"LoadOC_ListTraces()"
	PopupMenu Aligned_betapop,pos={40,59},size={223,21},bodyWidth=159,title="aligned beta:"
	PopupMenu Aligned_betapop,mode=1,popvalue="none",value= #"LoadOC_ListTraces()"
	PopupMenu Aligned_deltapop,pos={38,85},size={225,21},bodyWidth=159,title="aligned delta:"
	PopupMenu Aligned_deltapop,mode=1,popvalue="none",value= #"LoadOC_ListTraces()"
	Button MakeGraph,pos={54,165},size={159,30},proc=MakeGraph_but,title="Make Graph (for aligned OC)"
	Slider Adjust_slider,pos={4,111},size={260,52},proc=AdjustOC_slider
	Slider Adjust_slider,limits={-1,1,0},value= 0,vert= 0
	Button SaveOC_but,pos={129,213},size={135,36},proc=SaveOC_button,title="Save Optical Constants"
	SetVariable OCName,pos={4,222},size={119,16},title="Name"
	SetVariable OCName,value= root:OpticalConstants:Calcs:name
	setdatafolder foldersave
End
function sparcewaves(accuracy)
	variable accuracy
	string traces = tracenamelist("",";",1)
	wave xwave = xwavereffromtrace("",stringfromlist(0,traces))
	variable ntraces = itemsinlist(traces)
	make /free /n=(ntraces) /wave ywaves = tracenametowaveref("",stringfromlist(p,traces))
	make /wave /free /n=(ntraces) outputwaves, dwaves
	variable j, num=dimsize(ywaves,0), k
	string outputwavename, dwavename
	make /free /n=(num) lastlevel, direction=accuracy
	if(!waveexists(xwave))
		wave ywave = ywaves[0]
		duplicate /free ywave, xwave
		xwave=x
	endif
	differentiate xwave /d=tempwave
	variable minstep = wavemin(tempwave)
	
	make /n=3 /o xwaveout
	variable xmax = wavemax(xwave)
	variable xmin = wavemin(xwave)
	xwaveout[0] = xmin
	xwaveout[1] = .5*xmax + .5*xmin
	xwaveout[2] = xmax
	
	variable depth = ln((xmax-xmin)/minstep)/ln(2)
	
	for(j=0;j<num;j+=1)
		wave jwave = ywaves[j]
		outputwavename = uniquename(nameofwave(jwave),1,0)
		duplicate jwave, $outputwavename
		outputwaves[j] = $outputwavename
		appendtograph outputwaves[j] vs xwaveout
		modifygraph mode($outputwavename)=4
		
	endfor
	duplicate /free jwave, difwave, pwave
	// add endpoints and center point in xwaveout

	// itterate until maximum difference is below accuracy
	variable newpnts, xpnts, cnt=0
	do
		// calculate difference between original wave and interpolated new value
		newpnts = 0
		for(k=0;k<num;k+=1)
			wave ywave = ywaves[k]
			wave outwave = outputwaves[k]
			redimension /n=(dimsize(xwaveout,0)) outwave
			outwave = interp(xwaveout,xwave,ywave)
			difwave = abs(ywave - interp(xwave, xwaveout, outwave))
			findlevels /p /d=locations /q  difwave, accuracy
			if(dimsize(locations,0)==0)
				continue
			endif
			locations = xwave[locations]
			xpnts = dimsize(xwaveout,0)
			make /o /free /n=(xpnts,dimsize(locations,0)) addpoints
			// if there is a location between two data points that are more than one away from eachother, then bisect that section
			addpoints[0,xpnts-2][] = locations[q] >= xwaveout[p] && locations[q] < xwaveout[p+1]
			imagetransform sumAllRows, addpoints
			wave W_sumRows
			W_sumRows = W_sumRows? 1 : 0 // before which point should we add a point (in which gap)
				//  if w_sumrows[j] is 1, it means to add a point between j and j+1 where j runs from 0 to the second to last value in w_sumrows
			newpnts = sum(W_sumRows) // how many points to add
			for(j=xpnts-2; j>=0; j-=1) // go through the gaps backwards (starting at the last possible insert point, (xpnts-2) because the last point needs to remain
				if(W_sumRows[j]) // do we add a point to this gap? 
					insertpoints j+1 , 1 , xwaveout // add a point after point j
					xwaveout[j+1] = .5*xwaveout[j] + .5*xwaveout[j+2] // 
				endif
			endfor
		endfor
		for(k=0;k<num;k+=1)
			wave ywave = ywaves[k]
			wave outwave = outputwaves[k]
			redimension /n=(dimsize(xwaveout,0)) outwave
			outwave = interp(xwaveout,xwave,ywave)
		endfor
		doupdate
		cnt+=1
		if(cnt>depth+1)
			break
		endif
	while(newpnts) //points were addeed
	
	variable size = dimsize(xwaveout,0) 
	for(k=0;k<num;k+=1)
		wave ywave = ywaves[k]
		wave outwave = outputwaves[k]
		redimension /n=(size) outwave
		outwave = interp(xwaveout,xwave,ywave)
	endfor
	print size
	
end