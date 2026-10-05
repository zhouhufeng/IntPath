package musculus;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import tools.*;

public class main {
	
	public static String orgs = "musculus";
	
	//These HashMap are containers for the genes and gene pairs extraction results.
	private static HashMap<String,HashMap<String,String>> 
	cphge, cphgp, kphge, kphgp, wphge, wphgp;
		
	//These HashMap are containers for the genes and gene pairs normalization results.
	private static HashMap<String,HashMap<String,String>> 
	cge, cgp, kge, kgp, wge, wgp ;	
	
	private static HashMap<String,HashMap<String,Vector<String>>> 
	kgro,wgro,wprgro;
	
	private static HashMap<String,Vector<String>> idmp;
	
	public static void main(String[] args) throws Exception{
		
		Calendar c1 = Calendar.getInstance();
		
		extraction();
		
		normalization();
		
		statisticsBeforeInt();
		
		integration();
		
		statisticsAfterInt();
		
		preparePPIN();
		
		Calendar c2 = Calendar.getInstance();

		System.out.println( "Congratulations! All IntPath updating process in "+orgs+" is now successfully finished.\n"+
				"Running time of whole progrom takes："+(double)(c2.getTimeInMillis()-c1.getTimeInMillis())/1000+"s");
			
	}



	private static void preparePPIN() throws IOException{
		// This method is used to prepare the PPIN used in calculating the distances between pathways.
		
		String outmmu = orgs+File.separator+"source"+File.separator+"PPI"+File.separator+orgs+"STRING";		
		
		/*
		int thrldmmu = 750;
		
		conts.isoSTRINGPPI(orgs, "10090.", outmmu,thrldmmu);
		*/
		
		String cleanppi = orgs+File.separator+"normalized"+File.separator+"PPI"+File.separator+orgs+"STRING";	
		
		HashMap<String,String> ppimp = conts.PPImaping(idmp,outmmu,cleanppi);
		
	}




	private static void statisticsBeforeInt() {
		
		System.out.println("\nThe average number of genes and gene pairs per pathway before and after KEGG normalization.");
		/**/	
		calc kex = new calc();
		kex.AverGENnGPR(kphge, kphgp);		
		/**/	
		calc knm = new calc();		
		knm.AverGENnGPR(kge, kgp);
		
		System.out.println("\nThe average number of genes and gene pairs per pathway before and after BioCyc normalization.");
		/**/
		calc cex = new calc();
		cex.AverGENnGPR(cphge, cphgp);
		/**/
		calc cnm = new calc();		
		cnm.AverGENnGPR(cge, cgp);		
		
		System.out.println("\nThe average number of genes and gene pairs per pathway before and after WikiPathway normalization.");
		/**/		
		calc wex = new calc();		
		wex.AverGENnGPR(wphge, wphgp);
		/**/
		calc wnm = new calc();		
		wnm.AverGENnGPR(wge, wgp);
		
		
		//----Following is to calculate the overlap of genes and gene pairs among different databases;
		
		System.out.println("To calculate the overlap of genes between different databases.\n");
		
		calc.compDBGenes(kge,"KEGG",cge,"BioCyc");
		
		calc.compDBGenes(wge,"WikiPathways",cge,"BioCyc");
		
		calc.compDBGenes(kge,"KEGG",wge,"WikiPathways");
		
		
		System.out.println("To calculate the overlap of gene pairs between different databases.\n");
		
		calc.compDBGenePairs(kgp,"KEGG",cgp,"BioCyc");
		
		calc.compDBGenePairs(wgp,"WikiPathways",cgp,"BioCyc");
		
		calc.compDBGenePairs(kgp,"KEGG",wgp,"WikiPathways");
		
		
		//-----Calculate the specific pathway gene pairs and genes overlap------------
		
		//System.out.println("debug"+kge.size()+"\t"+kgp.size());
		
        HashMap<String,String> kpthge = kge.get("Citrate cycle (TCA cycle)");
        
        HashMap<String,String> kpthgp = kgp.get("Citrate cycle (TCA cycle)");
        
        System.out.println(
        		"number of genes in TCA pathway in KEGG"+kpthge.size()+
        		
        		"number of gene pairs in TCA pathway in KEGG"+kpthgp.size()
        		
        		);
        
        HashMap<String,String> wpthge = wge.get("TCA Cycle");
        
        HashMap<String,String> wpthgp = wgp.get("TCA Cycle");
        
        System.out.println(
        		"number of genes in TCA pathway in WikiPathway"+wpthge.size()+
        		
        		"number of gene pairs in TCA pathway in WikiPathway"+wpthgp.size()
        		
        		);
        
        HashMap<String,String> cpthge = cge.get("TCA cycle");
        
        HashMap<String,String> cpthgp = cgp.get("TCA cycle");
        
        System.out.println(
        		
        		"number of genes in TCA pathway in BioCyc"+cpthge.size()+
        		
        		"number of gene pairs in TCA pathway in BioCyc"+cpthgp.size()
        		
        		);
                
        
        System.out.println("To calculate the overlap of genes in TCA pathway between different databases.\n");
        
        calc.compGenes(kpthge,"KEGG TCA",wpthge,"Wiki TCA");
        
        calc.compGenes(kpthge,"KEGG TCA",cpthge,"Cyc TCA");
        
        calc.compGenes(cpthge,"Cyc TCA" ,wpthge,"Wiki TCA");


        //-------calculate the specific pathway gene pairs overlap ------
        System.out.println("To calculate the overlap of gene pairs in TCA pathway between different databases.\n");
        
        calc.compPairs(kpthgp,"KEGG TCA",wpthgp,"Wiki TCA");
        
        calc.compPairs(kpthgp,"KEGG TCA",cpthgp,"Cyc TCA");
        
        calc.compPairs(cpthgp,"Cyc TCA" ,wpthgp,"Wiki TCA");
		
	}

	private static void statisticsAfterInt() {		
				
		System.out.println("The average number of genes and gene pairs per pathway in the integrated pathways." );
		
		HashMap<String,HashMap<String,String>> onlyintgen = Integration.intGENOnly;
		
		HashMap<String,HashMap<String,String>> onlyintgpr = Integration.intGPROnly;
		
		System.out.println("The number of integrated pathways in pathway-genes:" +onlyintgen.size()+
				"\nThe number of integrated pathways in pathway-genepairs: "+onlyintgpr.size());
		
		
		
		System.out.println("\nThe average number of genes and gene pairs per pathway in the KEGG pathways after integration." );
		
		HashMap<String,HashMap<String,String>> rke = Integration.resKge;
		
		HashMap<String,HashMap<String,String>> rkp = Integration.resKgp;
		
		System.out.println("\nThe number of rest KEGG pathways in pathway-genes:" +rke.size()+
				"\nThe number of rest KEGG pathways in pathway-genepairs: "+rkp.size());
				
		
		
		System.out.println("\nThe average number of genes and gene pairs per pathway in the BioCyc pathways after integration." );
		
		HashMap<String,HashMap<String,String>> rce = Integration.resCge;// zero size
		
		HashMap<String,HashMap<String,String>> rcp = Integration.resCgp;
		
		System.out.println("\nThe number of rest BioCyc pathways in pathway-genes:" +rce.size()+
				"\nThe number of rest BioCyc pathways in pathway-genepairs: "+rcp.size());			
		
		
		
		System.out.println("\nThe average number of genes and gene pairs per pathway in the WikiPathways after integration." );
		
		HashMap<String,HashMap<String,String>> rwe = Integration.resWge;// zero size
		
		HashMap<String,HashMap<String,String>> rwp = Integration.resWgp;
		
		System.out.println("The number of rest pathways in WikiPathways pathway-genes:" +rwe.size()+
				"\nThe number of rest pathways in WikiPathways pathway-genepairs: "+rwp.size());
		
		
		System.out.println("\nThe average number of Integrated genes and gene pairs per pathway after Integration.");
		/**/	
		calc intall = new calc();
		intall.AverGENnGPR(onlyintgen, onlyintgpr);		
		
		/*	
		calc ain = new calc();		
		ain.AverGENnGPR(kge, kgp);
		*/		
		
		System.out.println("\nThe average number of KEGG genes and gene pairs per pathway after Integration.");
		/**/	
		calc kint = new calc();
		kint.AverGENnGPR(rke, rkp);
		
		/**/	
		calc kin = new calc();		
		kin.AverGENnGPR(kge, kgp);
		
		
		
		System.out.println("\nThe average number of BioCyc genes and gene pairs per pathway after Integration.");
		/**/
		calc cint = new calc();
		cint.AverGENnGPR(rce, rcp);		
		/**/
		calc cin = new calc();		
		cin.AverGENnGPR(cge, cgp);		
		
		
		System.out.println("\nThe average number of WikiPathways genes and gene pairs per pathway after Integration.");
		/**/		
		calc wint = new calc();		
		wint.AverGENnGPR(rwe, rwp);		
		/**/
		calc win = new calc();		
		win.AverGENnGPR(wge, wgp);
		
		
		//Calculate the number of pathways that is mutually contained in all three databases
		HashMap<String,Vector<String>> intpwyname = Integration.intpathnames;
		calc nameint = new calc();
		int thdbpwy = nameint.threeOverlapPWYs(intpwyname);
		System.out.println("The number of pathways that are mutually contained in three databases is: "+thdbpwy);
		
	}	

	private static void extraction() throws Exception{
		
		//-------------Extract BioCyc----------------------
		
		BioCyc.extract(orgs);
		
		//HashMap<String,HashMap<String,String>> 
		cphge = BioCyc.phge;
		
		//HashMap<String,HashMap<String,String>> 
		cphgp = BioCyc.phgp;
		
			
		calc cycgenes = new calc();
		double cphwygene = cycgenes.AverNUMperPWY(cphge);
		
		System.out.println("The number of pathways that have genes on average in Cyc: "+cphwygene);		
		
		calc cycgenpairs = new calc();
		double cphwygenpair = cycgenpairs.AverNUMperPWY(cphgp);
		
		System.out.println("The number of pathways that have gene pairs on average in Cyc: "+ cphwygenpair+"\n");

		/*	
		test extrBioCyc = new test(); 
		extrBioCyc.writeEveryThingOut(cphge);		
		extrBioCyc.writeEveryThingOut(cphgp);			
		*/

		
		//-------------Extract KEGG------------------------
				
		KEGG.extract(orgs);
		
		//HashMap<String,HashMap<String,String>> 
		kphge = KEGG.phge;
		
		//HashMap<String,HashMap<String,String>> 
		kphgp = KEGG.phgp;	
		
		kgro = KEGG.phgro;

		
		calc kegenes = new calc();
		double kphwygene = kegenes.AverNUMperPWY(kphge);
		System.out.println("The number of pathways that have genes on average in KEGG: "+kphwygene);	
		
		calc kegenpairs = new calc();
		double kphwygenpair = kegenpairs.AverNUMperPWY(kphgp);
		System.out.println("The number of pathways that have gene pairs on average in KEGG: "+kphwygenpair+"\n");		

		
		//This is just to double check if the work in extract KEGG is correct. 
		/*
		test extrKEGG = new test(); 
		extrKEGG.writeEveryThingOut(kphge);		
		extrKEGG.writeEveryThingOut(kphgp);  		
		*/ 

		
	
		//--------------Extract Wiki------------------------
		
		Wiki.extract(orgs);
		
		//HashMap<String,HashMap<String,String>> 
		wphge = Wiki.phge;
		
		//HashMap<String,HashMap<String,String>> 
		wphgp = Wiki.phgp;		
		
		//get group containers
		wprgro = Wiki.phgro;
		
		System.out.println("The num of WikiPathway gene Groups before normalization: "+ wprgro.size());
		
		calc wigenes = new calc();
		double wphwygene = wigenes.AverNUMperPWY(wphge);		
		System.out.println("The number of pathways that have genes on average in Mouse Wiki: "+wphwygene);	
		
		calc wigenpairs = new calc();
		double wphwygenpair = wigenpairs.AverNUMperPWY(wphgp);		
		System.out.println("The number of pathways that have gene pairs on average in Mouse Wiki: "+ wphwygenpair+"\n");				

		
		//This is just to double check if the work in extract KEGG is correct. 
		
		/*/*
		test extrWiki = new test(); 
		extrWiki.writeEveryThingOut(wphge);		
		extrWiki.writeEveryThingOut(wphgp);		
		*/		
		
	}

	private static void normalization() throws Exception{

		Normalize.buildmap(orgs);
		
		idmp = Normalize.mp;						
		
		createTables.ctGeneMapping(orgs,idmp);
				
		/*
		 * Write tow file
		 */
		String idmpingf = orgs+File.separator+"integrated"+File.separator+"Archive"+File.separator+"IDmappingFile";
		
		test testidmp = new test();
		
		testidmp.writempvToFile(idmp,idmpingf);
		

		
		//------------Normalize BioCyc genes and gene pairs-----------
		
		test excpairs = new test();		
		
		int excgp = excpairs.NumHashInHash(cphgp);
		
		int excge = excpairs.NumHashInHash(cphge);
		
		Normalize nc = new Normalize();
		
		//HashMap<String,HashMap<String,String>> 
		cge = nc.gemp(cphge, "BioCyc");		
		
		//HashMap<String,HashMap<String,String>> 
		cgp = nc.gpmp(cphgp);			
		
		test nomcpairs = new test();		
		
		int nomcgp = nomcpairs.NumHashInHash(cgp);
		
		int nomcge = nomcpairs.NumHashInHash(cge);
		
		
		
		System.out.println(
				"The number of BioCyc gene pairs before normalization: "+excgp+
				"\nThe number of BioCyc genes before normalization: "+excge+
				"\nThe number of BioCyc gene pairs after normalization: "+nomcgp+
				"\nThe number of BioCyc genes after normalization: "+nomcge
				);
		
		
		calc cycgenes = new calc();
		
		double cphgene = cycgenes.AverNUMperPWY(cge);
		
		System.out.println("The number of pathways that have genes on average in Nomarlized Cyc: "+cphgene);			
		
		calc cycgenpairs = new calc();
		
		double cphgenpair = cycgenpairs.AverNUMperPWY(cgp);
		
		System.out.println("The number of pathways that have gene pairs on average in Nomarlized Cyc: "+ cphgenpair+"\n");		
		
		//This is just to double check if the work in extract BioCyc is correct. 
		
		test normBioCyc = new test(); 
		/*
		normBioCyc.writeEveryThingOut(cge);		
		normBioCyc.writeEveryThingOut(cgp);		
		*/
		String ncge = orgs+File.separator+"normalized"+File.separator+"BioCyc"+File.separator+orgs+"BioCycNormPthGEN";				
		
		normBioCyc.writeGenesToFile(cge, ncge);
		
		String ncgp = orgs+File.separator+"normalized"+File.separator+"BioCyc"+File.separator+orgs+"BioCycNormPthGPR";
		
		normBioCyc.writeGenePairsToFile(cgp, ncgp, "BioCyc");
		
		
		//------------Normalize KEGG genes and gene pairs-----------		
		
		test exkpairs = new test();		
		
		int exkgp = exkpairs.NumHashInHash(kphgp);
		
		int exkge = exkpairs.NumHashInHash(kphge);
		
		
		Normalize nk = new Normalize();
		
		//HashMap<String,HashMap<String,String>> 
		kge = nk.gemp(kphge,"KEGG");
		
		//HashMap<String,HashMap<String,String>> 
		kgp = nk.gpmp(kphgp);	
		
		test nomkpairs = new test();		
		
		int nomkgp = nomkpairs.NumHashInHash(kgp);
		
		int nomkge = nomkpairs.NumHashInHash(kge);
		
		System.out.println(
				"The number of KEGG gene pairs before normalization: "+exkgp+
				"\nThe number of KEGG genes before normalization: "+exkge+
				"\nThe number of KEGG gene pairs after normalization: "+nomkgp+
				"\nThe number of KEGG genes after normalization: "+nomkge
				);		
		/**/		
		calc kegenes = new calc();
		double kphgene = kegenes.AverNUMperPWY(kge);		
		System.out.println("The number of pathways that have genes on average in Nomarlized  KEGG: "+kphgene);			
		
		calc kegenpairs = new calc();
		double kphgenpair = kegenpairs.AverNUMperPWY(kgp);		
		System.out.println("The number of pathways that have gene pairs on average in Nomarlized  KEGG: "+ kphgenpair+"\n");		
		
		//This is just to double check if the work in extract KEGG is correct. 		
				
		test normKEGG = new test(); 
		/*
		normKEGG.writeEveryThingOut(kge);		
		normKEGG.writeEveryThingOut(kgp);		
		*/
		String nkge = orgs+File.separator+"normalized"+File.separator+"KEGG"+File.separator+orgs+"KEGGNormPthGEN";				
		normKEGG.writeGenesToFile(kge, nkge);
		
		String nkgp = orgs+File.separator+"normalized"+File.separator+"KEGG"+File.separator+orgs+"KEGGNormPthGPR";
		normKEGG.writeGenePairsToFile(kgp, nkgp, "KEGG");
		
		
		//------------Normalize Wiki genes and gene pairs-----------
		
		test exwpairs = new test();		
		
		int exwgp = excpairs.NumHashInHash(wphgp);
		
		int exwge = excpairs.NumHashInHash(wphge);	
		
		Normalize nw = new Normalize();
		
		//HashMap<String,HashMap<String,String>> 
		wge = nw.gemp(wphge,"WikiPathways");
		
		//HashMap<String,HashMap<String,String>> 
		wgp = nw.gpmp(wphgp);
		
		wgro =nw.groupmp(wprgro);
		
		System.out.println("The num of WikiPathway gene Groups after normlization: "+ wgro.size());
		
		test nomwpairs = new test();		
		
		int nomwgp = nomcpairs.NumHashInHash(wgp);
		
		int nomwge = nomcpairs.NumHashInHash(wge);
		
		System.out.println(
				"The number of Wiki gene pairs before normalization: "+exwgp+
				"\nThe number of Wiki genes before normalization: "+exwge+
				"\nThe number of Wiki gene pairs after normalization: "+nomwgp+
				"\nThe number of Wiki genes after normalization: "+nomwge
				);			

		/**/
		calc wigenes = new calc();
		double wphgene = wigenes.AverNUMperPWY(wge);
		System.out.println(	"The number of pathways that have genes on average in  Wiki: "+wphgene);			
		
		calc wigenpairs = new calc();
		double wphgenpair = wigenpairs.AverNUMperPWY(wgp);		
		System.out.println("The number of pathways that have gene pairs on average in  Wiki: "+ wphgenpair+"\n");		
		
		//This is just to double check if the work in extract KEGG is correct. 
		
		test normWiki = new test(); 
		/*
		normWiki.writeEveryThingOut(wge);		
		normWiki.writeEveryThingOut(wgp);		
		*/
		String nwge = orgs+File.separator+"normalized"+File.separator+"WikiPathways"+File.separator+orgs+"WikiPathwaysNormPthGEN";				
		normWiki.writeGenesToFile(wge, nwge);
		
		String nwgp = orgs+File.separator+"normalized"+File.separator+"WikiPathways"+File.separator+orgs+"WikiPathwaysNormPthGPR";
		normWiki.writeGenePairsToFile(wgp, nwgp, "WikiPathways");
		
		
		/*
		 * The following codes are just testing if there are one to more mapping exists in the current mapping file. 
		 */
		
		HashMap<String,Vector<String>> tmp = nw.testmp;
		
		test tsidmp = new test();
		
		String testid = orgs+File.separator+"integrated"+File.separator+"Archive"+File.separator+orgs+"IDMappingDebugs";
		
		tsidmp.writempvToFile(tmp,testid);	
		
		/*
		 * The following steps are use to output the groups in pathway. 
		 * passing parameter, pw, source database tag, and group container.
		 * output into a file. 
		 */
		
		String groupf = orgs+File.separator+"integrated"+File.separator+"IntPathData"+File.separator+orgs+"GroupGenes";
		
		PrintWriter gppw = new PrintWriter(new BufferedWriter(new FileWriter(groupf)));
		
		test.OutputGroups(kgro,"KEGG",gppw);
		
		test.OutputGroups(wgro,"WikiPathways",gppw);
		
		gppw.close();
		
	}

	private static void integration() throws Exception{
		
		Integration.combine(cge, cgp, kge, kgp, wge, wgp,orgs);
		
		
		
	}	
}
