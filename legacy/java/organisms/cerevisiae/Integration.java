package cerevisiae;

import java.io.*;
import java.sql.SQLException;
import java.util.*;
import java.util.regex.*;
import java.util.Map.*;
import tools.*;

public class Integration{
	
	
	public static HashMap<String,HashMap<String,String>> intGENOnly = new HashMap<String,HashMap<String,String>>();
	
	public static HashMap<String,HashMap<String,String>> intGPROnly = new HashMap<String,HashMap<String,String>>();
	
	public static HashMap<String, Vector<String>> intpathnames  = new HashMap<String, Vector<String>>();

	
	// combine is one of the main functions of the Integration class.
	public static void combine(HashMap<String, HashMap<String, String>> cge,HashMap<String, HashMap<String, String>> cgp,
			HashMap<String, HashMap<String, String>> kge,HashMap<String, HashMap<String, String>> kgp,
			HashMap<String, HashMap<String, String>> wge,HashMap<String, HashMap<String, String>> wgp,
			String orgs) throws IOException, ClassNotFoundException, SQLException{
		
		String statsf     = orgs+File.separator+"integrated"+File.separator+"Stats"+File.separator+"IntegrationStatistics";
		
		PrintWriter stpw  = new PrintWriter(new BufferedWriter(new FileWriter(statsf))); 
		
		conts BCcon = new conts();
		Vector<String> cGEpthv = BCcon.keyToVec(cge);
		Vector<String> cGPpthv = BCcon.keyToVec(cgp);		
		System.out.println("Num. pathways from BioCyc (Gene)\t"+cGEpthv.size()+"\tGenePairs\t"+cGPpthv.size());
		
		conts KEcon = new conts();
		Vector<String> kGEpthv = KEcon.keyToVec(kge);
		Vector<String> kGPpthv = KEcon.keyToVec(kgp);
		System.out.println("Num. pathways from KEGG (Gene)\t"+kGEpthv.size()+"\tGenePairs\t"+kGPpthv.size());
		
		conts WIcon = new conts();
		Vector<String> wGEpthv = WIcon.keyToVec(wge);
		Vector<String> wGPpthv = WIcon.keyToVec(wgp);		
		System.out.println("Num. pathways from WikiPathways (Gene)\t"+wGEpthv.size()+"\tGenePairs\t"+wGPpthv.size());	
		
		/* This is writting the pathway names from each of the source databases into the 
		 * Mysql table "pathway_names" of corresponding organims. 
		 */
		
		createTables.ctPathwayNames(cGEpthv, kGEpthv, wGEpthv,orgs);
		
		/*
		 * The intPWYGroup is a HashMap, with integrated pathway name as key, 
		 * Original pathway name<k> and source database<v> stored in a HashMap as value.
		 * Since names in pathway-genes(GEpthv) is more than that in pathway-genepairs(GPpth). 
		 * So the groupPathwayNames have done specificly for the pathway-genes(GEpthv) and pathway-genepairs(GPpth).
		 */		
		
		System.out.println("\nThe following are the comibining steps for the pathway-genes statistics:\n");
		
		stpw.println("\nThe following are the comibining steps for the pathway-genes statistics:\n");
		
		String relPthsGef = orgs+File.separator+"integrated"+File.separator+"ReltedPathNames"+File.separator+"RelPthNamsGEN";		
		
		HashMap<String, Vector<String>> relPthsGE = groupPathwayNames(cGEpthv, kGEpthv, wGEpthv,relPthsGef,stpw);				
				
		
		/* This is writting the pathway names from each of the source databases into the 
		 * Mysql table "pathway_names" of corresponding organims. 
		 */
		
		intpathnames.putAll(relPthsGE);
		
		createTables.ctRelatedPathways(relPthsGE,orgs);
		
		/*	
		System.out.println("\nThe following are the comibining steps for the pathway-genepairs statistics:\n");
		
		stpw.println("\nThe following are the comibining steps for the pathway-genepairs statistics:\n");
		
		String relPthsGpf = orgs+File.separator+"integrated"+File.separator+"ReltedPathNames"+File.separator+"RelPthNamsGPR";
		
		HashMap<String, Vector<String>> relPthsGP = groupPathwayNames(cGPpthv, kGPpthv, wGPpthv,relPthsGpf,stpw);
		*/		
		
		
		// The following is creating a Data structure to host the final output of the integrated pathway-genes.
		HashMap<String,HashMap<String,String>> IntPathGEN = BuildIntPathGEN(relPthsGE,kge,wge,cge,stpw);
		
		String intpathgef = orgs+File.separator+"integrated"+File.separator+"IntPathData"+File.separator+orgs+"IntPathGenes";
		

		test testint = new test();
		
		// Write the integrated pathway-genes out into file from the IntPathGEN data structure.
		testint.writeGenesToFile(IntPathGEN, intpathgef);
		
		//insert all the data of pathway-genes from the data structure into the MySQL table "pathway_genes";
		createTables.ctPathwayGenes(IntPathGEN, orgs);
		
		// The following is creating a Data structure to host the final output of the integrated pathway-genepairs.
		HashMap<String,HashMap<String,String>> IntPathGPR = BuildIntPathGPR(relPthsGE,kgp,wgp,cgp,stpw);
		
		String intpathGPR = orgs+File.separator+"integrated"+File.separator+"IntPathData"+File.separator+orgs+"IntPathGenePairs";				
		
		// Write the integrated pathway-genepairs out into file from the IntPathGPR data structure.
		testint.writeGenePairsToFile(IntPathGPR, intpathGPR);
		
		//insert all the data of pathway-genes from the data structure into two MySQL table "gene_pairs" and " pathway_allpairinfo";
		createTables.ctPathwayGenePairs(IntPathGPR, orgs);
				
		
		stpw.close();
	}
	

	private static HashMap<String, HashMap<String, String>> BuildIntPathGPR(HashMap<String, Vector<String>> relPth,
			HashMap<String, HashMap<String, String>> kg,HashMap<String, HashMap<String, String>> wg,
			HashMap<String, HashMap<String, String>> cg, PrintWriter stpw) {				
		
        System.out.println(
                "Before pathway-genepairs integration the pathway number in each of three databases:\n K  W   C:  "+"\n"+
                kg.size()+"\t"+wg.size()+"\t"+cg.size()
                );
        
        stpw.println(
                "Before pathway-genepairs integration the pathway number in each of three databases:\n K  W   C:  "+"\n"+
                kg.size()+"\t"+wg.size()+"\t"+cg.size()
                );
        
        int ctrmv = 0; // count the number of removed pathways in each databases;
		
		HashMap<String,HashMap<String,String>> IntPathGPR = new HashMap<String,HashMap<String,String>>();
		
		Set<String> intpwyset = relPth.keySet();
		
		for(String intpwy : intpwyset ){
			
			Vector<String> pwyv = relPth.get(intpwy);
			
			
			HashMap<String,String> intgps = new HashMap<String,String>();
			
			
			
			for(int i = 0; i < pwyv.size(); i++){
				
				String pwydb = pwyv.get(i);
				
				String[] pdb = pwydb.split("\\+");
				
				//HashMap<String,HashMap<String,String>> ge = null;
				
				HashMap<String,String> epgp = null;
				
				String pwy = pdb[0];
				
				String sdb = "";
				
				if(pdb[1].equalsIgnoreCase("K")){
					
					if(kg.containsKey(pwy)){
						
						ctrmv++;
						
						epgp = kg.get(pwy);
						
					}										
					
					
					sdb = "KEGG";
					
					kg.remove(pwy);
					
					
					
				}else if (pdb[1].equalsIgnoreCase("C")){
				
					if(cg.containsKey(pwy)){
						
						ctrmv++;
						
						epgp = cg.get(pwy);
						
					}					
					
					sdb = "BioCyc";
					
					cg.remove(pwy);
																				
					
				}else if (pdb[1].equalsIgnoreCase("W")){
					
					if(wg.containsKey(pwy)){
						
						epgp = wg.get(pwy);
						
						ctrmv++;
						
					}
															
					sdb = "WikiPathways";
					
					wg.remove(pwy);										
					
				}else{
					
					System.out.println("Warning! In BuildIntPathGPR Some original pthway names do not have + tags");
					
				}					
				
				if(epgp!=null){
					
					//intgps.putAll(epgp);
					
					
					Set<String> epgpset = epgp.keySet();
					
					for(String gps: epgpset){
																					
						String intrel = "";
						
						String intdb  = "";
						
						String srel = epgp.get(gps);
						
						if(intgps.containsKey(gps)){

							
							String intreldb = intgps.get(gps);
							
							String[] reldb = intreldb.split("\t");
							
							intrel = reldb[0];
							
							for(int j = 1; reldb.length>1 && j<reldb.length;j++){
								
								intdb = intdb + " " + reldb[j];
								
							}
							
							if(srel.equalsIgnoreCase(reldb[0])){
								
								intdb = intdb + " " +sdb;
								
								intgps.put(gps, intrel.trim()+"\t"+intdb.trim());
								
								
							}else{
								
//								if(sdb.equalsIgnoreCase("KEGG")){
//									
//									intrel = srel;
//									
//								}
								if(!intrel.contains(srel)){
									
									intrel = intrel+" "+srel;
									
								} 
								
								intdb = intdb + " " +sdb;
								
								intgps.put(gps, intrel.trim()+"\t"+intdb.trim());
								
							}
							
//							if(!intsdb.contains(db)){
//								
//								intsdb = intsdb +"\t"+db;
//								
//								intgps.put(ges,intsdb);
//								
//							}
							
						}else{
							
							intgps.put(gps,srel.trim()+"\t"+sdb.trim());
							
						}
						
					}
					//*/
					
				}												
				
			}
			
			if(intgps.size()>0){
				
				IntPathGPR.put(intpwy,intgps);
				
			}
			
		}		     
		
		/*//The codes below is used to check the difference between my implement methods and putAll by Java.
		test tputall = new test();
		//java put all
		int tj = tputall.NumHashInHash(IntPathGPR);
		
		System.out.println("Checking....\nThe total number of gene pairs in IntPathGPR is: "+tj);
		*/
		
        intGPROnly.putAll(IntPathGPR);
		
        System.out.println(
                "After pathway-genepairs integration the pathway number in each of three databases:\n K  W   C G "+"\n"+
                kg.size()+"\t"+wg.size()+"\t"+cg.size()+"\t"+IntPathGPR.size()
                );
        
        stpw.println(
                "After pathway-genepairs integration the pathway number in each of three databases:\n K  W   C G "+"\n"+
                kg.size()+"\t"+wg.size()+"\t"+cg.size()+"\t"+IntPathGPR.size()+"\n"+
                "The number of removed pathways: "+ctrmv
                );
        
        System.out.println("The number of removed pathways: "+ctrmv+
        		"\nThen put all the pathway-genepairs remaining in three databases into the IntPathGEN...");
        
        //getRestGENmap(IntPathGEN,kg);
        
        //AddAllMP(IntPathGPR,kg,"KEGG");
        
        IntPathGPR = addAllMAP(IntPathGPR,kg,"KEGG");
        
        System.out.println("The number of pathways IntPath (pathway-genepairs) after get all rest of KEGG pathways: "+IntPathGPR.size()+". " +
        		"The number of pathways in KEGG didn't integrated before:"+kg.size());
        stpw.println("The number of pathways IntPath (pathway-genepairs) after get all rest of KEGG pathways: "+IntPathGPR.size()+". " +
        		"The number of pathways in KEGG didn't integrated before:"+kg.size());
        
        
        //AddAllMP(IntPathGPR,cg,"BioCyc");
        
        IntPathGPR = addAllMAP(IntPathGPR,cg,"BioCyc");
        
        System.out.println("The number of pathways IntPath (pathway-genepairs) after get all rest of BioCyc pathways: "+IntPathGPR.size()+". " +
        		"The number of pathways in BioCyc didn't integrated before:"+cg.size());
        stpw.println("The number of pathways IntPath (pathway-genepairs) after get all rest of BioCyc pathways: "+IntPathGPR.size()+". " +
        		"The number of pathways in BioCyc didn't integrated before:"+cg.size());
        
        
        //AddAllMP(IntPathGPR,wg,"WikiPathways");
        
        IntPathGPR = addAllMAP(IntPathGPR,wg,"WikiPathways");        
        
        System.out.println("The number of pathways IntPath (pathway-genepairs) after get all rest of WikiPathways: "+IntPathGPR.size()+". " +
        		"The number of pathways in KEGG didn't integrated before:"+wg.size());
        stpw.println("The number of pathways IntPath (pathway-genepairs) after get all rest of WikiPathways: "+IntPathGPR.size()+". " +
        		"The number of pathways in KEGG didn't integrated before:"+wg.size());
		
        
		test trt = new test();
		//java put all
		int rt = trt.NumHashInHash(IntPathGPR);		
		
		System.out.println("Checking..return and non-return..\n The total number of gene pairs in IntPathGPR is: "+rt);		
        
		return IntPathGPR;
	}


	private static HashMap<String, HashMap<String, String>> addAllMAP(HashMap<String, HashMap<String, String>> intPathGPR,
			HashMap<String, HashMap<String, String>> gp, String sdb) {
		
		Set<String> pgpset = gp.keySet();
		
		for(String pwy : pgpset){
			
			HashMap<String,String> intgp = new HashMap<String,String>();
			
			HashMap<String,String> gpmp  = gp.get(pwy);
			
			Set<String> pset = gpmp.keySet();
			
			for(String ps: pset){
				
				String rel = gpmp.get(ps);
				
				String reldb = rel+"\t"+sdb;
				
				intgp.put(ps, reldb);
				
			}
			
			intPathGPR.put(pwy, intgp);
			
			
		}
		
		return intPathGPR;
	}


	private static void AddAllMP(HashMap<String, HashMap<String, String>> intPathGPR,
			HashMap<String, HashMap<String, String>> gp, String sdb) {
		
		Set<String> pgpset = gp.keySet();
		
		for(String pwy : pgpset){
			
			HashMap<String,String> intgp = new HashMap<String,String>();
			
			HashMap<String,String> gpmp  = gp.get(pwy);
			
			Set<String> pset = gpmp.keySet();
			
			for(String ps: pset){
				
				String rel = gpmp.get(ps);
				
				String reldb = rel+"\t"+sdb;
				
				intgp.put(ps, reldb);
				
			}
			
			intPathGPR.put(pwy, intgp);
			
			
		}
		
	}


	private static HashMap<String, HashMap<String, String>> BuildIntPathGEN(HashMap<String, Vector<String>> relPth, 
			HashMap<String, HashMap<String, String>> kg, HashMap<String, HashMap<String, String>> wg, 
			HashMap<String, HashMap<String, String>> cg, PrintWriter stpw) {
		
        System.out.println(
                "Before pathway-genes integration the pathway number in each of three databases:\n K  W   C:  "+"\n"+
                kg.size()+"\t"+wg.size()+"\t"+cg.size()
                );
        
        stpw.println(
                "Before pathway-genes integration the pathway number in each of three databases:\n K  W   C:  "+"\n"+
                kg.size()+"\t"+wg.size()+"\t"+cg.size()
                );
        
        int ctrmv = 0; // count the number of removed pathways in each databases;
		
		HashMap<String,HashMap<String,String>> IntPathGEN = new HashMap<String,HashMap<String,String>>();
		
		Set<String> intpwyset = relPth.keySet();
		
		for(String intpwy : intpwyset ){
			
			Vector<String> pwyv = relPth.get(intpwy);
			
			
			HashMap<String,String> intgens = new HashMap<String,String>();
			
			
			
			for(int i = 0; i < pwyv.size(); i++){
				
				String pwydb = pwyv.get(i);
				
				String[] pdb = pwydb.split("\\+");
				
				//HashMap<String,HashMap<String,String>> ge = null;
							
				HashMap<String,String> epge = null;
				
				String pwy = pdb[0];
				
				if(pdb[1].equalsIgnoreCase("K")){
															
					epge = kg.get(pwy);
					
					kg.remove(pwy);
					
					ctrmv++;
					
				}else if (pdb[1].equalsIgnoreCase("C")){
				
					epge = cg.get(pwy);
					
					cg.remove(pwy);
					
					ctrmv++;
					
				}else if (pdb[1].equalsIgnoreCase("W")){
					
					epge = wg.get(pwy);
					
					wg.remove(pwy);
					
					ctrmv++;
					
				}else{
					
					System.out.println("Warning! In BuildIntPathGEN Some original pthway names do not have + tags");
					
				}					
				
				if(epge!=null){
					
					Set<String> epgeset = epge.keySet();
					
					for(String ges: epgeset){
						
						String db = epge.get(ges);
						
						if(intgens.containsKey(ges)){
							
							String intsdb = intgens.get(ges);
							
							if(!intsdb.contains(db)){
								
								intsdb = intsdb +" "+db;
								
								intgens.put(ges,intsdb);
								
							}
							
						}else{
							
							intgens.put(ges,db);
							
						}
						
					}
					
				}												
				
			}
			
			if(intgens.size()>0){
				
				IntPathGEN.put(intpwy,intgens);
				
			}
			
		}		       

        intGENOnly.putAll(IntPathGEN);
        
        System.out.println(
                "After pathway-genes integration the pathway number in each of three databases:\n K  W   C G "+"\n"+
                kg.size()+"\t"+wg.size()+"\t"+cg.size()+"\t"+IntPathGEN.size()
                );
        
        stpw.println(
                "After pathway-genes integration the pathway number in each of three databases:\n K  W   C G "+"\n"+
                kg.size()+"\t"+wg.size()+"\t"+cg.size()+"\t"+IntPathGEN.size()+"\n"+
                "The number of removed pathways: "+ctrmv
                );
        
        System.out.println("The number of removed pathways: "+ctrmv+
        		"\nThen put all the pathway-genes remaining in three databases into the IntPathGEN...");
        
        //getRestGENmap(IntPathGEN,kg);
        
        IntPathGEN.putAll(kg);
        
        System.out.println("The number of pathways IntPath (pathway-genes) after get all rest of KEGG pathways: "+IntPathGEN.size()+". " +
        		"The number of pathways in KEGG didn't integrated before:"+kg.size());
        stpw.println("The number of pathways IntPath (pathway-genes) after get all rest of KEGG pathways: "+IntPathGEN.size()+". " +
        		"The number of pathways in KEGG didn't integrated before:"+kg.size());
        
        
        IntPathGEN.putAll(cg);
        
        System.out.println("The number of pathways IntPath (pathway-genes) after get all rest of BioCyc pathways: "+IntPathGEN.size()+". " +
        		"The number of pathways in BioCyc didn't integrated before:"+cg.size());
        stpw.println("The number of pathways IntPath (pathway-genes) after get all rest of BioCyc pathways: "+IntPathGEN.size()+". " +
        		"The number of pathways in BioCyc didn't integrated before:"+cg.size());
        
        
        IntPathGEN.putAll(wg);
        
        System.out.println("The number of pathways IntPath (pathway-genes) after get all rest of WikiPathways: "+IntPathGEN.size()+". " +
        		"The number of pathways in KEGG didn't integrated before:"+wg.size());
        stpw.println("The number of pathways IntPath (pathway-genes) after get all rest of WikiPathways: "+IntPathGEN.size()+". " +
        		"The number of pathways in KEGG didn't integrated before:"+wg.size());
		
		return IntPathGEN;
	}

	/*//didin't implement because output the number increases each time can also trace if there is any replace of mappings happened. 

	private static void getRestGENmap( HashMap<String, HashMap<String, String>> intPathGEN,
			HashMap<String, HashMap<String, String>> kg) {
		
		
		
	}
	*/

	private static HashMap<String, Vector<String>> groupPathwayNames(
			Vector<String> cpthv, Vector<String> kpthv,Vector<String> wpthv, 
			String relPthsGe,PrintWriter stpw)throws IOException {
		
		/*
		 * Bug log(2011-Nov.7 11:30): 
		 * comparePathwayNames have problem failed to identify the related name, 
		 * The way to process integrated pathway name is wrong.
		 *  
		 */
		
		HashMap<String, HashMap<String, String>> pthgroup = new HashMap<String, HashMap<String, String>>();
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(relPthsGe)));
		
		pw.println("The pair-wise related pathway pairs: ");
		
		Vector<String[]> kwv = comparePathwayNames  (kpthv, "K",  wpthv, "W", pw);
		
		Vector<String[]> kcv = comparePathwayNames  (kpthv, "K",  cpthv, "C", pw);
		
		Vector<String[]> cwv = comparePathwayNames  (cpthv, "C",  wpthv, "W", pw);
		
		Vector<String[]> kkv = compareInPathwayNames(kpthv, "K", pw);
		
		Vector<String[]> wwv = compareInPathwayNames(wpthv, "W", pw);
		
		Vector<String[]> ccv = compareInPathwayNames(cpthv, "C", pw);			
		
//		Vector<String[]> kkv = comparePathwayNames(kpthv, "K",  kpthv, "K", pw);
//		
//		Vector<String[]> wwv = comparePathwayNames(wpthv, "W",  wpthv, "W", pw);
//		
//		Vector<String[]> ccv = comparePathwayNames(cpthv, "C",  cpthv, "C", pw);				
		
        System.out.println("Pairwise Results of pathways that can be merged: \n" +
        		" KEGG-Wiki\tKEGG-BioCyc\twiki-BioCyc\tKEGG-KEGG\tBioCyc-BioCyc \t wiki-wiki \n"
                +kwv.size()+" \t "+kcv.size()+" \t "+cwv.size()+" \t "+kkv.size()+" \t "+ccv.size()+" \t "+wwv.size());
		
        stpw.println("Pairwise Results of pathways that can be merged: \n" +
        		" KEGG-Wiki\tKEGG-BioCyc\twiki-BioCyc\tKEGG-KEGG\tBioCyc-BioCyc \t wiki-wiki \n"
                +kwv.size()+" \t "+kcv.size()+" \t "+cwv.size()+" \t "+kkv.size()+" \t "+ccv.size()+" \t "+wwv.size());
		
        HashMap<String,String> dsjtmp = new HashMap<String,String>();
        
        
        buildDisjointSet(kwv, dsjtmp);
        
        //System.out.println("Pass in kwv (size:"+kwv.size()+")the num of nods disjoint set: "+dsjtmp.size());
        
        buildDisjointSet(kcv, dsjtmp);
        
        //System.out.println("Pass in kcv (size:"+kcv.size()+")the num of nods disjoint set: "+dsjtmp.size());
        
        buildDisjointSet(cwv, dsjtmp);
        
        //System.out.println("Pass in cwv (size:"+cwv.size()+")the num of nods disjoint set: "+dsjtmp.size());
        
        buildDisjointSet(kkv, dsjtmp);
        
        //System.out.println("Pass in kkv (size:"+kkv.size()+")the num of nods disjoint set: "+dsjtmp.size());
        
        buildDisjointSet(wwv, dsjtmp);
        
        //System.out.println("Pass in wwv (size:"+wwv.size()+")the num of nods disjoint set: "+dsjtmp.size());
        
        buildDisjointSet(ccv, dsjtmp);
        
        //System.out.println("Pass in ccv (size:"+ccv.size()+") the num of nods disjoint set: "+dsjtmp.size());
        
        System.out.println("No. of pathways from 3 source pathway Data Bases are meged:"+dsjtmp.size());
        
        stpw.println("No. of pathways from 3 source pathway Data Bases are meged:"+dsjtmp.size());
        
                
        HashMap<String, Vector<String>> initgroup = getPthGroups(dsjtmp);
        
        HashMap<String, Vector<String>> pathgroup = identifyIntPWYname(initgroup);
         
        
        pw.println("The integrated pathway name and correspond original pathway name: ");
                
        test.writeRelPWYgroup(pathgroup, pw);
        
        System.out.println("The number of integrated pathway groups: "+pathgroup.size()+"\t Double check\t"+ initgroup.size());
        
        stpw.println("The number of integrated pathway groups: "+pathgroup.size()+"\t Double check\t"+ initgroup.size());
                
        
        
        pw.close();
        
		return pathgroup;
	}

    private static HashMap<String, Vector<String>> identifyIntPWYname(
			HashMap<String, Vector<String>> initgroup) {
    	
    	HashMap<String, Vector<String>> pathgroup = new HashMap<String, Vector<String>>();
    	
    	Set parentset = initgroup.keySet();
    	
    	for(Object parentO : parentset){
    		
    		String parent = parentO.toString();
    		
    		Vector nodsv  = initgroup.get(parent);
    		
			//set an initial long value here, while replace by shorter strings.    			
			String upathdb = "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" +
					"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" +
					"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";
    		
    		for(int n =0; n<nodsv.size();n++){
    			
    			String nod = (String) nodsv.get(n);
    			
    			if(nod.length()<upathdb.length()){
    				
    				upathdb = nod;
    				
    			}
    			
    		}
    		
    		String[] pdb = upathdb.split("\\+"); 
    		//pdb[0] shortest pathway name, pdb[1] source of database where shortest name come from"K,C,W"
    		
    		String pthOriName = pdb[0];
        	
        	String un = processIntPathNames(pthOriName);
    		
    		//String un = pthOriName;
    		
            if(pathgroup.containsKey(un)){
            	
            	System.out.println("Alert! Something wrong here, two integrated pathway name in different pathway groups.");
            	
            }else{
            	
            	pathgroup.put(un, nodsv);
            	
            }
    		    		    	
    		
    	}
    	
		return pathgroup;
	}


//	private static String processIntPathNames(String pthOriName) {
//		
//		//System.out.println(pthOriName);
//		
//		String un = "";
//		
//    	String[] pnas = pthOriName.split(" ");
//		
//    	
//        for(int x= 0;x<pnas.length;x++){
//        	
//        	Boolean tell = ckNameParts(pnas[x]);
//        	
//        	if (tell){
//        		
//        		un = un + pnas[x]+" ";
//        		
//        	}
//        	
//        }
//        		
//
//		return un;
//	}


	private static HashMap<String, Vector<String>> getPthGroups(HashMap<String, String> dsjtmp) {
		
    	HashMap<String, Vector<String>> mp = new HashMap<String, Vector<String>>();
    	
    	Set nodset = dsjtmp.keySet();    	    	
    	
    	for(Object nodO:nodset){
    		
    		String nod = nodO.toString();
    		
    		String parent = recursivelygetParent(nod, dsjtmp);
    		
    		Vector<String> v = new Vector<String>();
    		
    		if(mp.containsKey(parent)){
    			
    			v = mp.get(parent);
    			
    			if(!v.contains(nod)){
    				
    				v.add(nod);
    				
    				mp.put(parent, v);
    				
    			}
    			
    		}else{
    			
				v.add(nod);
				
				mp.put(parent, v);
    			
    		}
    		
    	}
    	
		return mp;
	}


	private static void buildDisjointSet(Vector<String[]> kwv,HashMap<String,String> dsjtmp) {

        for(int x = 0; x<kwv.size();x++){
            
        	String[] kw = kwv.get(x);
            
            String pas = kw[0]+"+"+kw[1];
            
            String pbs = kw[2]+"+"+kw[3];
            
            if(!dsjtmp.containsKey(pas)&&!dsjtmp.containsKey(pbs)){
            	
            	dsjtmp.put(pas,pas);
                
            	dsjtmp.put(pbs,pas);
            	
            }else if (dsjtmp.containsKey(pas) && !dsjtmp.containsKey(pbs)){
                
            	dsjtmp.put(pbs,pas);
            	
            }else if (!dsjtmp.containsKey(pas)&& dsjtmp.containsKey(pbs)){
                
            	dsjtmp.put(pas,dsjtmp.get(pbs));
                //System.out.println("Be careful, this step might have errors.It is all right");
            }else{
            	
                String ap =  recursivelygetParent(pas,dsjtmp);
                
                String bp =  recursivelygetParent(pbs,dsjtmp);
                
                dsjtmp.put(ap, bp);
            }
        }
    	
	}


    private static String recursivelygetParent(String nod,HashMap<String,String> dsjtmp) {
        if(nod.equalsIgnoreCase(dsjtmp.get(nod))){
            return nod;
        }else{
            return recursivelygetParent(dsjtmp.get(nod),dsjtmp);
        }
    }


	private static Vector<String[]> comparePathwayNames (Vector<String> kv, String kdb,
    		Vector<String> wv, String wdb,PrintWriter pw) throws IOException {
    	
    	Vector<String[]> kwv = new Vector<String[]>();

        for (int k = 0; k < kv.size(); k++) {

            //double[] kwmax = {0.0, 0.0};
            
            double threld = 0.0;
            
            double aligSc = 0.0; 
            
            String[] kws = {"a", "a", "a", "a"};
            
            String pki = kv.get(k);
            
            String pk  = pki.trim();

            for (int w = 0; w < wv.size(); w++) {
            	
            	String wki = wv.get(w);
            	
            	String wk  = wki.trim();
            	
//                if (kv.get(k) == null || wv.get(w) == null) {
//                    continue;
//                }
                
                sequenceAlignment sa = new sequenceAlignment(pk, wk);
                
                if (sa.thres > threld) {//previously there was the bug. 
                	aligSc = sa.alsc;
                    threld = sa.thres;
                    kws[0] = pk;
                    kws[1] = kdb;
                    kws[2] = wk;
                    kws[3] = wdb;
                    
                    //System.out.println(kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
                }
            }

            //System.out.println( kws[1]+"\t"+kws[0]+"\t"+kwmax[0]+"\t"+kwmax[1]);

            if (aligSc > kws[0].length() - 1 || aligSc > kws[2].length() - 1 || threld > 0.91) {
                
            	//if(threld<0.5) continue;
            	boolean flag = ckmismatch(kws);	
            	
            	if (threld >= 0.5 && flag) {
            		
                    kwv.add(kws);
                    
                    pw.println(kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
                    
                    //System.out.println(flag+"\t"+kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
                }
            }

        }
        
        return kwv;
    }
	
    private static Vector<String[]> compareInPathwayNames(Vector<String> kv,
			String kdb, PrintWriter pw) throws IOException {
    	
    	Vector<String[]> kwv = new Vector<String[]>();

        for (int k = 0; k < kv.size(); k++) {

            //double[] kwmax = {0.0, 0.0};
            
            double threld = 0.0;
            
            double aligSc = 0.0; 
            
            String[] kws = {"a", "a", "a", "a"};
            
            String pki = kv.get(k);
            
            String pk  = pki.trim();

            for (int w = k+1; w < kv.size(); w++) {
            	
            	String wki = kv.get(w);
            	
            	String wk  = wki.trim();
            	
//                if (kv.get(k) == null || wv.get(w) == null) {
//                    continue;
//                }
                
                sequenceAlignment sa = new sequenceAlignment(pk, wk);
                
                if (sa.thres > threld) {//previously there was the bug. 
                	aligSc = sa.alsc;
                    threld = sa.thres;
                    kws[0] = pk;
                    kws[1] = kdb;
                    kws[2] = wk;
                    kws[3] = kdb;
                    
                    //System.out.println(kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
                }
            }

            //System.out.println( kws[1]+"\t"+kws[0]+"\t"+kwmax[0]+"\t"+kwmax[1]);

            if (aligSc > kws[0].length() - 1 || aligSc > kws[2].length() - 1 || threld > 0.91) {
                
            	//if(threld<0.5) continue;
            	
            	boolean flag = ckmismatch(kws);
            	
            	if (threld >= 0.5 && flag) {
                    
            		kwv.add(kws);
                    
            		pw.println(kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
            		
            		//System.out.println(flag+"\t"+kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
            		
                }
            	
            	
            }

        }
        
        return kwv;
	}
	


	private static String processIntPathNames(String pthOriName) {
		
		//System.out.println(pthOriName);
		
		String un = "";
		
    	String[] pnas = pthOriName.split(" ");
		
    	
        for(int x= 0;x<pnas.length;x++){
        	
        	Boolean tell = ckNameParts(pnas[x]);
        	//System.out.println("un:"+un);
        	
        	if (tell){
        		
        		
        		pnas[x] = SpecialNameReplace(pnas[x]);
        		
        		//System.out.println(un);
        		
        		un = un + pnas[x]+" ";
        		
        	}
        	
        }
        		

		return un;
	}


	private static Boolean ckNameParts(String parts) {
		
		String[] dlst = {
				"I","","II","V","III","IV","VI","Type","X"
				//,"(",")","even","number","branch","ADP-D-Glucose","pathway"
				,"adenosine","keratan","heparan","-","globo","Small","",""//This line is human specific
				,"saturated",",","","","","","",""
				};
		
		for (int i = 0; i < dlst.length;i++){
			
			if(parts.equalsIgnoreCase(dlst[i])) return false;
			
			
		}
		
		
		return true;
	}


	private static String SpecialNameReplace(String un) {
		
		String nun = "";
		
		if (un.equalsIgnoreCase("O-Glycan")){
			
			nun = "Glycan";
			
			return nun;
			
		}else if(un.equalsIgnoreCase("IL-2")){
			
			nun ="Interleukin";						
			
			return nun;
			
		}else if(un.equalsIgnoreCase("ubiquinone-8")){//
			
			nun ="ubiquinone";						
			
			return nun;
			
		}else if(un.equalsIgnoreCase("palmitate")){
			
			nun ="palmitoleate and palmitate";						
			
			return nun;
			
		}else if(un.equalsIgnoreCase("glycolysis")){
			
			nun ="Glycolysis and Gluconeogenesis";						
			
			return nun;
			
		}else if(un.equalsIgnoreCase("Deoxyribose")){
			
			nun ="Ribose and Deoxyribose";						
			
			return nun;
			
		}else if(un.equalsIgnoreCase("leucine")){
			
			nun ="Valine, leucine and isoleucine";						
			
			return nun;
			
		}else if(un.equalsIgnoreCase("spermine")){
			
			nun ="spermine, spermidine";						
			
			return nun;
			
		}else{
			
			return un;
			
		}
		
	}
	
	private static boolean ckmismatch(String[] kws) {
		
//		boolean flag = ckmismatch(kws);	
		
		String[][] err ={ 
			{"NOD","Toll"},
			{"Linoleic","Lipoic"},
			{"T cell","B cell"},
			{"EPO","TOR"},
			{"L-cysteine","lysine"},
			{"spermine","serine"},
			{"serotonin","serine"},
			{"Steroid","thyroid"},
//			{"isoleucine","leucine"},
			{"dermatan","heparan"},
			{"ribonucleotides","deoxyribonucleotides"},
			{"guanosine","adenosine"},
			{"sulfation","oxidation"},
			{"Lysine","Glycine"},
			{"Serine","homoserine"},
			{"galactose","Lactose"},
//			{"siroheme","Heme"},
			{"alanine biosynthesis","Valine Biosynthesis"},
//			{"Leucine","Isoleucine"},
			{"valine biosynthesis","alanine biosynthesis"},
			{"oleate","folate"},
			{"glutaredoxin","thioredoxin"},	
			{"homoserine","homocysteine"},
			{"valine degradation","alanine degradation"},
//			{"methionine","threonine"},
			{"phenylalanine biosynthesis","alanine biosynthesis"},
			{"glutathione","glutamine"},
			{"S-adenosylmethionine","Methionine"},
			{"Glycerolipid","Glycerophospholipid"},
			{"guanine, xanthine","adenine, hypoxanthine"},
//			{"",""},
		};
		
		
		for(int i = 0; i<err.length;i++){
			
			if(kws[0].contains(err[i][0])&&kws[2].contains(err[i][1])) {
	
				System.out.println("Mis-matches of related pathway names handled!\n"+
						kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3]);
				
				return false;
								
			}else if(kws[2].contains(err[i][0])&&kws[0].contains(err[i][1])){
				
				System.out.println("Mis-matches of related pathway names handled!\n"+
						kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3]);
				
				return false;
			}
		
		}		
		
		return true;
		
	}	
 
    
}
