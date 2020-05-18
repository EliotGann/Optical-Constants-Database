#pragma rtGlobals=1	// Use modern global access method.
Menu "Macros"
	submenu "Optical Constants database"
		"Upload NEXAFS to Optical Constant DB", /Q, Execute/P/Q "OpenLoadOCPanel()"
			help={"Load calculated Delta and Betas from a graph into a database of Optical Constants"}
		"Load Optical Constants DB", /Q, Execute/P/Q "LoadOC()"
			help={"Loads opticalconstants.txt from User Procedures"}
	end
End
function loadOC_old()
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
			saveallOCs()
			//Save/O/G/W/U={0,0,1,0}/B /p=opticalconstants wavelist("*",";","") as "OpticalConstants.txt"
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


function testHDF5OCwrite(path, filename, EN_WAVE, isodelta,isobeta[paradelta,parabeta,perpdelta,perpbeta])
	string path, filename
	wave en_wave, isodelta, isobeta,paradelta,parabeta,perpdelta,perpbeta
	variable unaligned = paramIsDefault(paradelta)
	variable fileID
	HDF5CreateFile/O /P=$path fileID  as fileName
		if(v_flag)
		print "Failed to write hdf file"
		return -1
	endif
	// get metadata from the wave itself
	string facility = stringbykey("facility",note(en_wave))
	string instrument = stringbykey("instrument",note(en_wave))
	string technique = stringbykey("technique",note(en_wave))
	string technique_details = stringbykey("technique_details",note(en_wave))
	string normalization_method = stringbykey("normalization_method",note(en_wave))
	string publication_doi = stringbykey("publication_doi",note(en_wave))
	string submitter = stringbykey("submitter",note(en_wave))
	string submitter_email = stringbykey("submitter_email",note(en_wave))
	string corresponding_author = stringbykey("corresponding_author",note(en_wave))
	string corresponding_author_email = stringbykey("corresponding_author_email",note(en_wave))
	string chemical_formula = stringbykey("chemical_formula",note(en_wave))
	string samplename = stringbykey("samplename",note(en_wave))
	string description = stringbykey("description",note(en_wave))
	variable density = numberbykey("density",note(en_wave))
	variable thickness = numberbykey("thickness",note(en_wave))
	
	if(density*thickness*0 != 0) // wave does not contain metadata, go for manual entry
		dfref foldersave = getdataFolderDFR()
		setdatafolder root:
		newdatafolder /o/s opticalconstants
		// check for default values
		svar /z gdescription, gfacility, ginstrument,gtechnique, gtechnique_details, gnormalization_method, gpublication_doi
		svar /z gsubmitter, gsubmitter_email, gcorresponding_author, gcorresponding_author_email, gchemical_formula
		nvar /z gthickness,gdensity
		if(!svar_exists(gfacility))
			string /g gfacility = "NSLS-II"
		endif
		if(!svar_exists(ginstrument))
			string /g ginstrument = "SST-1 RSoXS"
		endif
		if(!svar_exists(gtechnique))
			string /g gtechnique = "TEY"
		endif
		if(!svar_exists(gtechnique_details))
			string /g gtechnique_details = "measured with a channeltron detector with 1900V acceleration bias"
		endif
		if(!svar_exists(gnormalization_method))
			string /g gnormalization_method = note(isodelta)
		endif
		if(!svar_exists(gpublication_doi))
			string /g gpublication_doi = "N/A"
		endif
		if(!svar_exists(gsubmitter))
			string /g gsubmitter = "Eliot Gann"
		endif
		if(!svar_exists(gsubmitter_email))
			string /g gsubmitter_email = "eliot.gann@nist.gov"
		endif
		if(!svar_exists(gcorresponding_author_email))
			string /g gcorresponding_author_email = "eliot.gann@nist.gov"
		endif
		if(!svar_exists(gcorresponding_author))
			string /g gcorresponding_author = "Eliot Gann"
		endif
		if(!svar_exists(gchemical_formula))
			string /g gchemical_formula = "C8H8"
		endif
		string /g gsamplename = replacestring(".oc",filename,"")
		if(!svar_exists(gdescription))
			string /g gdescription = note(isodelta)
		endif
		if(!nvar_exists(gdensity))
			variable /g gdensity = 1.2
		endif
		if(!nvar_exists(gthickness))
			variable /g gthickness = 100
		endif
		setdatafolder foldersave
		// set default values to the input values
		facility = gfacility
		instrument = ginstrument
		technique = gtechnique
		technique_details = gtechnique_details
		normalization_method = gnormalization_method
		publication_doi = gpublication_doi
		submitter = gsubmitter
		submitter_email = gsubmitter_email
		corresponding_author = gcorresponding_author
		corresponding_author_email = gcorresponding_author_email
		chemical_formula = gchemical_formula
		samplename = gsamplename
		description = gdescription
		density = gdensity
		thickness = gthickness	
		
		// prompt for information for sample
		Prompt samplename, "Sample name (simple name): "		// Set prompt for x param
		Prompt description, "Sample Description (as detailed as you want): "		// Set prompt for y param
		Prompt chemical_formula, "Chemical Formual (eg C1H4): "		// Set prompt for y param
		Prompt density, "Enter density of sample (g/cm^3): "		// Set prompt for x param
		Prompt thickness, "thickness of sample (nm): "		// Set prompt for y param
		Prompt submitter, "Submitter Name (you): "		// Set prompt for y param
		Prompt submitter_email, "Submitter Email address: "		// Set prompt for y param
		
		doprompt	"Sample Information" samplename, description,chemical_formula, density, thickness, submitter, submitter_email
		
		//prompt for information about measurement
		Prompt facility, "Facility where data was collected: "		// Set prompt for x param
		Prompt instrument, "Instrument used to collect data: "		// Set prompt for y param
		Prompt technique, "NEXAFS Technique (TEY/PEY/etc): "		// Set prompt for x param
		Prompt technique_details, "NEXAFS details: "		// Set prompt for x param
		Prompt normalization_method, "Describe the normalization method: "		// Set prompt for y param
		Prompt publication_doi, "Publication DOI:  "		// Set prompt for x param
		Prompt corresponding_author, "Corresponding Author:  "		// Set prompt for x param
		Prompt corresponding_author_email, "Corresponding Author Email address: "		// Set prompt for y param
		
		doprompt "Measurement/Origin Information" facility,instrument,technique,technique_details,normalization_method, publication_doi,corresponding_author, corresponding_author_email
		
		// save these values as default values
		gfacility = facility
		ginstrument = instrument
		gtechnique = technique
		gtechnique_details = technique_details
		gnormalization_method = normalization_method
		gsubmitter = submitter
		gsubmitter_email = submitter_email
		gpublication_doi = publication_doi
		gcorresponding_author = corresponding_author
		gcorresponding_author_email = corresponding_author_email
		gchemical_formula = chemical_formula
		gsamplename = samplename
		gdescription = description
		gdensity = density
		gthickness = thickness	
	endif
	// add root attributes
	HDFstringAttr(fileid,"default","entry","/")
	HDFstringAttr(fileid,"file_name","filename","/")
	HDFvarAttr(fileid,"file_time",datetime,"/")
	HDFstringAttr(fileid,"creator","QANT","/")
	HDFstringAttr(fileid,"creator_version","1.11","/")
	HDFstringAttr(fileid,"NeXus_version","4.3.0","/")
	HDFvarAttr(fileid,"HDF5_version",HDF5LibraryInfo(0),"/")
	HDFvarAttr(fileid,"HFD5XOP_version",HDF5LibraryInfo(0),"/")
	//add entry group
	variable NXEntryID
	HDF5CreateGroup fileID , "entry" , NXEntryID
	HDFstringAttr(fileID,"NX_class","NXentry","entry")
	HDFstringAttr(fileID,"default","dielectric_function","entry")
	HDFstringDS(NXEntryID,"title",description)
	//add data group
	variable nxdataid
	HDF5CreateGroup NXEntryID , "dielectric_function" , nxdataid
	
	HDFstringAttr(NXEntryID,"NX_class","NXdata","dielectric_function")
	HDFstringAttr(NXEntryID,"signal","delta,beta","dielectric_function")
	HDFstringAttr(NXEntryID,"delta_axes","energy,orientation","dielectric_function")
	HDFstringAttr(NXEntryID,"beta_axes","energy,orientation","dielectric_function")
	
	
	duplicate /free en_wave,energy
	HDF5SaveData /IGOR=0/O /Z energy, nxdataid, "energy"
	HDFstringAttr(nxdataid,"units","electronvolt","energy")
	HDFstringAttr(nxdataid,"long_name","Energy [Ev]","energy")
	
	if(!unaligned)
		//HDFstringAttr(NXEntryID,"orientation_type","uniaxial","dielectric_function")
		make /free orientations = {{0,0,0},{1,0,0},{0,1,1}}
		make /free /t orientation_type = {"isotropic","face-on parallel","face-on perpendicular"}
		make /free /t orientation_details = {"isotropic","direction normal to plane of aromatic core","average of directions planar to aromatic core"}
		//matrixtranspose orientations
		concatenate /free {isodelta, paradelta, perpdelta}, deltawave
		concatenate /free {isobeta, parabeta, perpbeta}, betawave
	else
		//HDFstringAttr(NXEntryID,"orientation_type","isotropic","dielectric_function")
		make /free orientations = {0,0,0}
		make /free /t orientation_type = {"isotropic"}
		make /free /t orientation_details = {"isotropic"}
		duplicate /free isodelta, deltawave
		duplicate /free isobeta, betawave
	endif
	
	HDF5SaveData /IGOR=0/O /Z orientations, nxdataid, "orientation_vector"
	HDF5SaveData /IGOR=0/O /Z /A="orientation_type" orientation_type, nxdataid, "orientation_vector"
	HDF5SaveData /IGOR=0/O /Z /A="orientation_details" orientation_details, nxdataid, "orientation_vector"
	HDF5SaveData /IGOR=0/O /Z deltawave, nxdataid, "delta"
	HDFstringAttr(nxdataid,"units","","delta")
	HDFstringAttr(nxdataid,"long_name","delta","delta")
	HDF5SaveData /IGOR=0/O /Z betawave, nxdataid, "beta"
	HDFstringAttr(nxdataid,"units","","beta")
	HDFstringAttr(nxdataid,"long_name","delta","beta")
	
	//add origin group
	variable originid
	HDF5CreateGroup NXEntryID , "optical_constants_origin" , originID
	HDFstringDS(originID,"type","NEXAFS")
	HDFvarAttr(NXEntryID,"time",datetime,"optical_constants_origin")
	HDFstringAttr(originID,"facility",facility,"type")
	HDFstringAttr(originID,"instrument",instrument,"type")
	HDFstringAttr(originID,"technique",technique,"type")
	HDFstringAttr(originID,"technique_details",technique_details,"type")
	HDFstringAttr(originID,"normalization_method",normalization_method,"type")
	HDFstringAttr(NXEntryID,"publication_doi",publication_doi,"optical_constants_origin")
	HDFstringAttr(NXEntryID,"corresponding_author",corresponding_author,"optical_constants_origin")
	HDFstringAttr(NXEntryID,"corresponding_author_email",corresponding_author_email,"optical_constants_origin")
	HDFstringAttr(NXEntryID,"submitter",submitter,"optical_constants_origin")
	HDFstringAttr(NXEntryID,"submitter_email",submitter_email,"optical_constants_origin")
	HDFstringAttr(NXEntryID,"database_name","NIST Optical Constants DB","dielectric_function")
	HDFstringAttr(NXEntryID,"database_accept_hash","untracked","dielectric_function")
	
	
	variable sampleid
	HDF5CreateGroup NXEntryID , "sample" , sampleid
	HDFstringAttr(NXEntryID,"NX_class","NXsample","sample")
	HDFstringDS(sampleid,"type","sample")
	HDFstringDS(sampleid,"name",samplename)
	HDFstringDS(sampleid,"description",description)
	HDFstringDS(sampleid,"chemical_formula",chemical_formula)
	HDFvarDS(sampleid,"density",density)
	HDFvarDS(sampleid,"thickness",thickness)
	HDFstringAttr(sampleid,"units","g/cm^3","density")
	HDFstringAttr(sampleid,"units","nm","thickness")
	
	HDF5CloseFile /z fileid
	

end
function HDFstringAttr(variable fileid, string attrname,string attr,string loc)
	// Write an attribute to the dataset
	if(strlen(loc)<1)
		loc = "/"
	endif
	Make /FREE /T /N=1 datasetAttribute = attr
	HDF5SaveData/IGOR=0 /A=attrname datasetAttribute, fileID, loc
end
function HDFstringDS(variable locid, string dsname,string sval)
	// Write a string dataset
	Make /FREE /T /N=1 dataset = sval
	HDF5SaveData/IGOR=0 dataset, locid, dsname
end
function HDFvarAttr(variable fileid, string attrname,variable attr,string loc)
	// Write an attribute to the dataset
	if(strlen(loc)<1)
		loc = "/"
	endif
	Make /FREE /N=1 datasetAttribute = attr
	HDF5SaveData/IGOR=0 /A=attrname datasetAttribute, fileID, loc
end
function HDFvarDS(variable locid, string dsname,variable val)
	// Write a variable dataset
	Make /FREE /N=1 dataset = val
	HDF5SaveData/IGOR=0 dataset, locid, dsname
end
function HDFwaveAttr(variable fileid, string attrname,wave attr,string loc)
	// Write an attribute to the dataset
	if(strlen(loc)<1)
		loc = "/"
	endif
	HDF5SaveData/IGOR=0 /A=attrname attr, fileID, loc
end


function saveallOCs()
	string matlist = listOCs()
	variable i, dosave
	string mat
	for(i=0;i<itemsinlist(matlist);i++)
		mat = stringfromlist(i,matlist)
		if(OCsexist(mat)==2)
			doalert /T="Save?" 1, "Aligned Material "+mat+" found.  save?"
			if(v_flag==1)
				wave enw = $("root:OpticalConstants:energy_"+mat)
				wave dw = $("root:OpticalConstants:delta_"+mat)
				wave bw = $("root:OpticalConstants:beta_"+mat)
				wave dwe = $("root:OpticalConstants:energy_"+mat+"perp")
				wave dwa = $("root:OpticalConstants:delta_"+mat+"perp")
				wave bwe = $("root:OpticalConstants:beta_"+mat+"para")
				wave bwa = $("root:OpticalConstants:beta_"+mat+"para")
				testhdf5OCwrite("opticalconstants",mat+".oc",enw,dw,bw,paradelta = dwa,parabeta = bwa,perpdelta = dwe,perpbeta= bwe)
			endif
		elseif(OCsExist(mat)==1)
			doalert /T="Save?" 1, "unaligned Material "+mat+" found.  save?"
			if(v_flag==1)
				wave enw = $("root:OpticalConstants:energy_"+mat)
				wave dw = $("root:OpticalConstants:delta_"+mat)
				wave bw = $("root:OpticalConstants:beta_"+mat)
				testhdf5OCwrite("opticalconstants",mat+".oc",enw,dw,bw)
			endif
		endif
	endfor
end

function loadOC()
	string foldersave = getdatafolder(1)
	setdatafolder root:
	newdatafolder /s/o opticalconstants
	//newpath /m="Select location of optical constant database: opticalconstants.txt " /O /Q opticalconstants
	NewPath/C/Z/Q/O opticalconstants SpecialDirPath("Igor Pro User files",0,0,0)+"User Procedures:Optical Constants"
	//GetFileFolderInfo /P=opticalconstants /Q /Z "opticalconstants.txt"
	string listoffiles = indexedFile(opticalconstants,-1,".oc")
	variable i
	string filename
	for(i=0;i<itemsinlist(listoffiles);i++)
		loadOCfile(filename)
	endfor
	setdatafolder foldersave
end

function loadOCfile(string filename)
	variable fileid
	
	HDF5OpenFile /R /z /P=opticalconstants fileid as filename
	if(v_flag)
		return -1
	endif
	newdatafolder /free loading
	
	variable NXentryID
	HDF5openGroup fileid,"entry",NXentryID
	HDF5LoadData NXentryID,"Title"
	variable nxdataid
	HDF5openGroup NXentryID,"dielectric_function",nxdataid
	HDF5loadData nxdataid, "energy"
	
	if(!unaligned)
		//HDFstringAttr(NXEntryID,"orientation_type","uniaxial","dielectric_function")
		make /free orientations = {{0,0,0},{1,0,0},{0,1,1}}
		make /free /t orientation_type = {"isotropic","face-on parallel","face-on perpendicular"}
		make /free /t orientation_details = {"isotropic","direction normal to plane of aromatic core","average of directions planar to aromatic core"}
		//matrixtranspose orientations
		concatenate /free {isodelta, paradelta, perpdelta}, deltawave
		concatenate /free {isobeta, parabeta, perpbeta}, betawave
	else
		//HDFstringAttr(NXEntryID,"orientation_type","isotropic","dielectric_function")
		make /free orientations = {0,0,0}
		make /free /t orientation_type = {"isotropic"}
		make /free /t orientation_details = {"isotropic"}
		duplicate /free isodelta, deltawave
		duplicate /free isobeta, betawave
	endif
	
	HDF5SaveData /IGOR=0/O /Z orientations, nxdataid, "orientation_vector"
	HDF5SaveData /IGOR=0/O /Z /A="orientation_type" orientation_type, nxdataid, "orientation_vector"
	HDF5SaveData /IGOR=0/O /Z /A="orientation_details" orientation_details, nxdataid, "orientation_vector"
	HDF5SaveData /IGOR=0/O /Z deltawave, nxdataid, "delta"
	HDFstringAttr(nxdataid,"units","","delta")
	HDFstringAttr(nxdataid,"long_name","delta","delta")
	HDF5SaveData /IGOR=0/O /Z betawave, nxdataid, "beta"
	HDFstringAttr(nxdataid,"units","","beta")
	HDFstringAttr(nxdataid,"long_name","delta","beta")
	
	//add origin group
	variable originid
	HDF5CreateGroup NXEntryID , "optical_constants_origin" , originID
	HDFstringDS(originID,"type","NEXAFS")
	HDFvarAttr(NXEntryID,"time",datetime,"optical_constants_origin")
	HDFstringAttr(originID,"facility",facility,"type")
	HDFstringAttr(originID,"instrument",instrument,"type")
	HDFstringAttr(originID,"technique",technique,"type")
	HDFstringAttr(originID,"technique_details",technique_details,"type")
	HDFstringAttr(originID,"normalization_method",normalization_method,"type")
	HDFstringAttr(NXEntryID,"publication_doi",publication_doi,"optical_constants_origin")
	HDFstringAttr(NXEntryID,"corresponding_author",corresponding_author,"optical_constants_origin")
	HDFstringAttr(NXEntryID,"corresponding_author_email",corresponding_author_email,"optical_constants_origin")
	HDFstringAttr(NXEntryID,"submitter",submitter,"optical_constants_origin")
	HDFstringAttr(NXEntryID,"submitter_email",submitter_email,"optical_constants_origin")
	HDFstringAttr(NXEntryID,"database_name","NIST Optical Constants DB","dielectric_function")
	HDFstringAttr(NXEntryID,"database_accept_hash","untracked","dielectric_function")
	
	
	variable sampleid
	HDF5CreateGroup NXEntryID , "sample" , sampleid
	HDFstringAttr(NXEntryID,"NX_class","NXsample","sample")
	HDFstringDS(sampleid,"type","sample")
	HDFstringDS(sampleid,"name",samplename)
	HDFstringDS(sampleid,"description",description)
	HDFstringDS(sampleid,"chemical_formula",chemical_formula)
	HDFvarDS(sampleid,"density",density)
	HDFvarDS(sampleid,"thickness",thickness)
	HDFstringAttr(sampleid,"units","g/cm^3","density")
	HDFstringAttr(sampleid,"units","nm","thickness")
	
	HDF5CloseFile /z fileid
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
end