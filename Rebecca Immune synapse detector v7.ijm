// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
//PV 2021Dec03
//Rebecca Immune Synapse detector
// --------------------
// --------------------
//Variables

TestImage01="D:/Image Analysis Workflows/Rebecca/CancerCells2.tif";
TestImage02="D:/Image Analysis Workflows/Rebecca/ImmuneCells2.tif";

CancerCellBGSub=500;
CancerCellKuwahara=11;
CancerCellGaussian=5;
CancerCellInteration=3;
CancerCellPixelCount=3;
CancerCell=0;
CancerCellMinimumThreshold =120;

ImmuneCellBGSub=5;
ImmuneCellKuwahara=5;
ImmuneCellGaussian=3;
ImmuneCellInteration=3;
ImmuneCellPixelCount=3;
ImmuneCellMinimumSize=200;
ImmuneCell=0;
ImmuneCellMinimumThreshold = 20;

MinimumOverlapArea=100;

OverlapingTimepoints=1;
OverlapingArea=100;
OverlapingAreaDilation=10;


// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
//Code starts below
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------

// ------------------------------------------------------------
//std initialise macro. also gets macro start time for reference
run("Collect Garbage");
run("Close All");		
run("Clear Results");		
roiManager("reset");		
run("Set Measurements...", "area redirect=None decimal=3");
if (isOpen("Log")) { selectWindow("Log"); run("Close");};		
if (isOpen("Summary")) { selectWindow("Summary"); run("Close");};
getDateAndTime(StartYear, StartMonth, StartDayOfWeek, StartDayOfMonth, StartHour, StartMinute, StartSecond, StartMsec);
setTool("zoom");
run("Options...", "iterations=1 count=1 black do=Nothing");
//exit;

// ------------------------------------------------------------
//edit analysis parameters

Dialog.create("Analysis parameters");		
Dialog.addMessage("Edit the analysis parameters below\n");	

Dialog.addMessage("\n");	
Dialog.addMessage("Cancer cells Parameters");	
Dialog.addNumber ("Background subtraction value", CancerCellBGSub );
Dialog.addNumber ("Smoothing (edge preserving)", CancerCellKuwahara );
Dialog.addNumber ("Bluring", CancerCellGaussian );
Dialog.addNumber ("Minimum Threshold Value", CancerCellMinimumThreshold );
Dialog.addNumber ("Close (Pixel)", CancerCellPixelCount );
Dialog.addNumber ("Close (Iteration)", CancerCellInteration );


Dialog.addMessage("\n");	
Dialog.addMessage("Immune cells Parameters");	
Dialog.addNumber ("Background subtraction value", ImmuneCellBGSub );
Dialog.addNumber ("Smoothing (edge preserving)", ImmuneCellKuwahara );
Dialog.addNumber ("Bluring", ImmuneCellGaussian );
Dialog.addNumber ("Minimum Threshold Value", ImmuneCellMinimumThreshold );
Dialog.addNumber ("Close (Pixel)", ImmuneCellPixelCount );
Dialog.addNumber ("Close (Iteration)", ImmuneCellInteration );
Dialog.addNumber ("Minimum Immune cell size", ImmuneCellMinimumSize );

Dialog.addMessage("\n");	
Dialog.addMessage("Other Parameters");	
Dialog.addNumber ("Minimum overlap area between Immune cell can cancer cell to be considered synapse", MinimumOverlapArea );
Dialog.addNumber ("Minimum overlap area between time points to be considered as synapse", OverlapingArea );
Dialog.addNumber ("Dilation of cells to detect continuity between timepoints", OverlapingAreaDilation );

Dialog.show();

CancerCellBGSub = Dialog.getNumber();
CancerCellKuwahara = Dialog.getNumber();
CancerCellGaussian = Dialog.getNumber();
CancerCellMinimumThreshold = Dialog.getNumber();
CancerCellInteration = Dialog.getNumber();
CancerCellPixelCount = Dialog.getNumber();

ImmuneCellBGSub = Dialog.getNumber();
ImmuneCellKuwahara = Dialog.getNumber();
ImmuneCellGaussian = Dialog.getNumber();
ImmuneCellMinimumThreshold = Dialog.getNumber();
ImmuneCellInteration = Dialog.getNumber();
ImmuneCellPixelCount = Dialog.getNumber();

MinimumOverlapArea = Dialog.getNumber();
OverlapingArea= Dialog.getNumber();
OverlapingAreaDilation=Dialog.getNumber();


//exit

// --------------------
//open files and rename channels

FileOpen= File.openDialog("Select File for analysis"); 	
OriginalName=substring(FileOpen, lastIndexOf(FileOpen, "\\")+1, lengthOf(FileOpen));
OriginalCoreName=substring(FileOpen, lastIndexOf(FileOpen, "\\")+1, lastIndexOf(FileOpen, "."));
OriginalExtnName=substring(FileOpen, lastIndexOf(FileOpen, "."), lengthOf(FileOpen));
OriginalDir=substring(FileOpen, 0, lastIndexOf(FileOpen, "\\")+1);
OutputFileDirectory=OriginalDir+OriginalCoreName+"_Output\\";	
File.makeDirectory(OutputFileDirectory);

//Open file using bioformats
run("Bio-Formats Importer", "open=FileOpen color_mode=Default rois_import=[ROI manager] view=[Standard ImageJ] stack_order=Default split_channels");					

//TestImage01= File.openDialog("Select Cancer Cells image File for analysis"); 	
//TestImage02= File.openDialog("Select Immune Cells image File for analysis"); 	
//open(TestImage01);		
//open(TestImage02);		

selectWindow(OriginalName+" - C=0");
rename("CancerCells");	run("Fire");

selectWindow(OriginalName+" - C=1");
rename("ImmuneCells");	run("Fire");

close(OriginalName+" - C=2");

//exit

// --------------------
//Immune cell processing
selectWindow("ImmuneCells");
run("Subtract Background...", "rolling=ImmuneCellBGSub sliding stack");
run("Kuwahara Filter", "sampling=ImmuneCellKuwahara stack");
run("Gaussian Blur...", "sigma=ImmuneCellGaussian stack");
setThreshold(ImmuneCellMinimumThreshold , 65535);
run("Convert to Mask", "method=Triangle background=Dark black");
run("Options...", "iterations=ImmuneCellInteration count=ImmuneCellPixelCount black do=Close stack");
run("Analyze Particles...", "size=ImmuneCellMinimumSize-Infinity pixel show=Masks stack");
run("Invert LUT");
rename("ImmuneCellsFiltered");


// --------------------
//Cancer cell processing
selectWindow("CancerCells");
run("Subtract Background...", "rolling=CancerCellBGSub sliding stack");
run("Kuwahara Filter", "sampling=CancerCellKuwahara stack");
run("Gaussian Blur...", "sigma=CancerCellGaussian stack");
setThreshold(CancerCellMinimumThreshold , 65535);
run("Convert to Mask", "method=Triangle background=Dark black");
run("Options...", "iterations=CancerCellInteration count=CancerCellPixelCount black do=Close stack");


// --------------------
//Synapse detection
imageCalculator("AND create stack", "CancerCells","ImmuneCellsFiltered");
rename("OverlapArea");
run("Analyze Particles...", "size=MinimumOverlapArea-Infinity pixel show=Masks stack");
run("Invert LUT");
rename("OverlapAreaSizeFiltered");

// --------------------
//Morpho reconstruction loop
selectWindow("ImmuneCellsFiltered");
getDimensions(ImmuneImageWidth, ImmuneImageHeight, ImmuneImageChannels, ImmuneImageSlices, ImmuneImageFrames);
run("Duplicate...", "title=ImmuneCellsWithSynapses duplicate");


for (i=1; i<=ImmuneImageSlices; i++)

	{
		selectWindow("ImmuneCellsFiltered");
		run("Duplicate...", "duplicate range=i-i use");
		rename("Temp_Original");


		selectWindow("OverlapAreaSizeFiltered");
		run("Duplicate...", "duplicate range=i-i use");
		rename("Temp_Overlap");
		
		run("Morphological Reconstruction", "marker=Temp_Overlap mask=Temp_Original type=[By Dilation] connectivity=8"); 

		rename("Temp_ImmuneCellsWithSynapses");
		run("Select All");
		run("Copy");

		selectWindow("ImmuneCellsWithSynapses");
		setSlice(i);
		run("Paste");

		close("Temp_Original");
		close("Temp_Overlap");
		close("Temp_ImmuneCellsWithSynapses");
	}


selectWindow("ImmuneCellsWithSynapses");
getDimensions(ImmuneImageWidth, ImmuneImageHeight, ImmuneImageChannels, ImmuneImageSlices, ImmuneImageFrames);
run("Duplicate...", "title=ImmuneCellsWithSustainedSynapsesForward duplicate");
run("Duplicate...", "title=ImmuneCellsWithSustainedSynapsesReverse duplicate");

//subroutine to detect persistant synapse as defined by cells with synapses overlapping with a certain area in teh next time point
//forward
for (i=1; i<=(ImmuneImageSlices-OverlapingTimepoints); i++)

	{
		//select first slice in test
		selectWindow("ImmuneCellsWithSynapses");
		setSlice(i);		
		run("Duplicate...", "title=temp1");

		//dilate to ctach nearby cells too. defined by variable
		run("Morphological Filters", "operation=Dilation element=Disk radius=OverlapingAreaDilation");
		rename("temp2");

		//select 2nd slice in test
		selectWindow("ImmuneCellsWithSynapses");
		setSlice(i+1);		
		run("Duplicate...", "title=temp3");

		//dilate to ctach nearby cells too. defined by variable
		run("Morphological Filters", "operation=Dilation element=Disk radius=OverlapingAreaDilation");
		rename("temp4");

		//look for overlap
		imageCalculator("AND create", "temp2","temp4");
		rename("temp5");

		//filter small overlaps out . defined in variable
		run("Analyze Particles...", "size=100-Infinity pixel show=Masks");
		run("Invert LUT");
		rename("temp6");

		//find original dilatd cell with overlap
		run("Morphological Reconstruction", "marker=temp6 mask=temp2 type=[By Dilation] connectivity=4");
		rename("temp7");

		//find original cell
		run("Morphological Reconstruction", "marker=temp7 mask=temp1 type=[By Dilation] connectivity=4");
		rename("temp8");		run("Select All");		run("Copy");

		//paste over original
		selectWindow("ImmuneCellsWithSustainedSynapsesForward");
		setSlice(i);		run("Paste");

		//find original dilatd cell with overlap
		run("Morphological Reconstruction", "marker=temp6 mask=temp4 type=[By Dilation] connectivity=4");
		rename("temp9");

		//find original cel
		run("Morphological Reconstruction", "marker=temp9 mask=temp3 type=[By Dilation] connectivity=4");
		rename("tempX");		run("Select All");		run("Copy");

		//paste over original	
		selectWindow("ImmuneCellsWithSustainedSynapsesForward");
		setSlice(i+1);		run("Paste");
		
		//cleanup
		close("temp1");close("temp2");close("temp3");close("temp4");close("temp5");
		close("temp6");close("temp7");close("temp8");close("temp9");close("tempX");
		
	}

//reverse

//reverse the order to chekc form other side
selectWindow("ImmuneCellsWithSynapses");
run("Reverse");

for (i=1; i<=(ImmuneImageSlices-OverlapingTimepoints); i++)

	{
		//select first slice in test
		selectWindow("ImmuneCellsWithSynapses");
		setSlice(i);		
		run("Duplicate...", "title=temp1");

		//dilate to ctach nearby cells too. defined by variable
		run("Morphological Filters", "operation=Dilation element=Disk radius=OverlapingAreaDilation");
		rename("temp2");

		//select 2nd slice in test
		selectWindow("ImmuneCellsWithSynapses");
		setSlice(i+1);		
		run("Duplicate...", "title=temp3");

		//dilate to ctach nearby cells too. defined by variable
		run("Morphological Filters", "operation=Dilation element=Disk radius=OverlapingAreaDilation");
		rename("temp4");

		//look for overlap
		imageCalculator("AND create", "temp2","temp4");
		rename("temp5");

		//filter small overlaps out . defined in variable
		run("Analyze Particles...", "size=100-Infinity pixel show=Masks");
		run("Invert LUT");
		rename("temp6");

		//find original dilatd cell with overlap
		run("Morphological Reconstruction", "marker=temp6 mask=temp2 type=[By Dilation] connectivity=4");
		rename("temp7");

		//find original cell
		run("Morphological Reconstruction", "marker=temp7 mask=temp1 type=[By Dilation] connectivity=4");
		rename("temp8");		run("Select All");		run("Copy");

		//paste over original
		selectWindow("ImmuneCellsWithSustainedSynapsesReverse");
		setSlice(i);		run("Paste");

		//find original dilatd cell with overlap
		run("Morphological Reconstruction", "marker=temp6 mask=temp4 type=[By Dilation] connectivity=4");
		rename("temp9");

		//find original cel
		run("Morphological Reconstruction", "marker=temp9 mask=temp3 type=[By Dilation] connectivity=4");
		rename("tempX");		run("Select All");		run("Copy");

		//paste over original	
		selectWindow("ImmuneCellsWithSustainedSynapsesReverse");
		setSlice(i+1);		run("Paste");
		
		//cleanup
		close("temp1");close("temp2");close("temp3");close("temp4");close("temp5");
		close("temp6");close("temp7");close("temp8");close("temp9");close("tempX");
		
	}

//reverse orders back to normal
selectWindow("ImmuneCellsWithSustainedSynapsesReverse");
run("Reverse");

selectWindow("ImmuneCellsWithSynapses");
run("Reverse");

//OR function to fill in cells in both directions
imageCalculator("OR create stack", "ImmuneCellsWithSustainedSynapsesForward","ImmuneCellsWithSustainedSynapsesReverse");
rename("ImmuneCellsWithSustainedSynapses");

//cleanup
close("ImmuneCellsWithSustainedSynapsesForward");close("ImmuneCellsWithSustainedSynapsesReverse");


//close("OverlapAreaSizeFiltered");

imageCalculator("Subtract create stack", "ImmuneCellsFiltered","ImmuneCellsWithSustainedSynapses");
rename("ImmuneCellsWithOutSynapses");

run("Merge Channels...", "c1=CancerCells c2=ImmuneCellsWithSustainedSynapses c3=ImmuneCellsWithOutSynapses keep");
RunVar=OutputFileDirectory+OriginalCoreName+"_Red_CancerCells_Green_ImmunecellsWithSustainedSynapse_Blue_ImmuneCellsWithOutSynapses";
saveAs("Tiff", RunVar);

close("CancerCells");
close("ImmuneCells");
close("ImmuneCellsWithSynapses");
//close("ImmuneCellsFiltered");
close("OverlapArea");
close("OverlapAreaSizeFiltered");

selectWindow("ImmuneCellsFiltered");
RunVar=OutputFileDirectory+OriginalCoreName+"_ImmuneCellsFiltered";
saveAs("Tiff", RunVar);

selectWindow("ImmuneCellsWithSustainedSynapses");
RunVar=OutputFileDirectory+OriginalCoreName+"_ImmuneCellsWithSustainedSynapses";
saveAs("Tiff", RunVar);

selectWindow("ImmuneCellsWithOutSynapses");
RunVar=OutputFileDirectory+OriginalCoreName+"_ImmuneCellsWithOutSynapses";
saveAs("Tiff", RunVar);



