/* Title: Debugging
 * Function: Java Program to Debugging 
 * 
 * Author: Hufeng Zhou
 * Time: July 3rd 2021
 * Version: v2
 */


package sapiens;
import java.util.HashMap;
import java.util.Set;
import java.util.Vector;
import java.io.*;
import java.util.*;
import tools.conts;
import tools.sequenceAlignment;
import tools.test;
/*
 * as we noticed the distinct differences between the old and new output files.
 * This bebuging program is used to tesing the integration process.  
 */
public class Debugging {
	
	public static void main(String args[]){
		
//		String[][] err ={ 
//				{"NOD","Toll"},
//				{"Linoleic","Lipoic"},
//				{"T cell","B cell"},
//				{"EPO","TOR"},
//				{"L-cysteine","lysine"},
//				{"spermine","serine"},
//				{"serotonin","serine"},
//			};
//		
//		for(int i = 0; i<err.length;i++){
//			for(int j = 0; j< err.length;j++){
//				
//			}
//			
//			System.out.println(err[i][0]+"\t"+err[i][1]);
//			
//		}
		
		String[] kws = {"NOD-like receptor signaling pathway","","Toll-like receptor signaling pathway",""};
		
		boolean flag = ckmismatch(kws);				
		
		if(flag){
			
			System.out.println("Good!");
			
			
		}else{
			
			System.out.println("mismatch!");
			
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
//			{"",""},
//			{"",""},
//			{"",""},
//			{"",""},
//			{"",""},
//			{"",""},
		};
		
		
		for(int i = 0; i<err.length;i++){
			
			if(kws[0].contains(err[i][0])&&kws[2].contains(err[i][1])) {
	
				System.out.println("Miss matches of related pathway names handled!\n"+
						kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3]);
				
				return false;
								
			}else if(kws[2].contains(err[i][0])&&kws[0].contains(err[i][1])){
				
				System.out.println("Miss matches of related pathway names handled!\n"+
						kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3]);
				
				return false;
			}
		
		}		
		
		return true;
		
	}	
	
	//public static PrintWriter bugpw = new PrintWriter(new BufferedWriter(new FileWriter("sapiens"+File.separator+"integrated"+File.separator+"Archive"+File.separator+"sapiensDebugs")));
	// combine is one of the main functions of the Integration class.
	public static void combine(HashMap<String, HashMap<String, String>> cge,HashMap<String, HashMap<String, String>> cgp,
			HashMap<String, HashMap<String, String>> kge,HashMap<String, HashMap<String, String>> kgp,
			HashMap<String, HashMap<String, String>> wge,HashMap<String, HashMap<String, String>> wgp,
			String orgs) throws IOException{
		
		PrintWriter bugpw = new PrintWriter(new BufferedWriter(new FileWriter("sapiens"+File.separator+"integrated"+File.separator+"Archive"+File.separator+"sapiensDebugs")));
		
		//Calendar c1 = Calendar.getInstance();
		
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
		
		/*
		 * The intPWYGroup is a HashMap, with integrated pathway name as key, 
		 * Original pathway name<k> and source database<v> stored in a HashMap as value.
		 * Since names in pathway-genes(GEpthv) is more than that in pathway-genepairs(GPpth). 
		 * So the groupPathwayNames have done specificly for the pathway-genes(GEpthv) and pathway-genepairs(GPpth).
		 */		
		
		System.out.println("\nThe following are the comibining steps for the pathway-genes statistics:\n");
		
		String relPthsGef = orgs+File.separator+"integrated"+File.separator+"ReltedPathNames"+File.separator+"RelPthNamsGEN";		
		
		HashMap<String, Vector<String>> relPthsGE = groupPathwayNames(cGEpthv, kGEpthv, wGEpthv,relPthsGef,stpw, bugpw);						
		
		
		
		System.out.println("\nThe following are the comibining steps for the pathway-genepairs statistics:\n");
		
		String relPthsGpf = orgs+File.separator+"integrated"+File.separator+"ReltedPathNames"+File.separator+"RelPthNamsGPR";
		
		HashMap<String, Vector<String>> relPthsGP = groupPathwayNames(cGPpthv, kGPpthv, wGPpthv,relPthsGpf,stpw, bugpw);
					
		
		String intpathGEN = orgs+File.separator+"integrated"+File.separator+"IntPathData"+File.separator+orgs+"IntPathGenes";
		
		String intpathGPR = orgs+File.separator+"integrated"+File.separator+"IntPathData"+File.separator+orgs+"IntPathGenePairs";
				
		
		stpw.close();
		
		bugpw.close();
		
//		Calendar c2 = Calendar.getInstance();
//
//		System.out.println("Running time it takes："+(c2.getTimeInMillis()-c1.getTimeInMillis())+"ms");
	}
	

	private static HashMap<String, Vector<String>> groupPathwayNames(
			Vector<String> cpthv, Vector<String> kpthv,Vector<String> wpthv, 
			String relPthsGe,PrintWriter stpw, PrintWriter bugpw)throws IOException {
		
		/*
		 * Bug log(2011-Nov.7 11:30): 
		 * comparePathwayNames have problem failed to identify the related name, 
		 * The way to process integrated pathway name is wrong.
		 *  
		 */
		
		HashMap<String, HashMap<String, String>> pthgroup = new HashMap<String, HashMap<String, String>>();
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(relPthsGe)));
				
		pw.println("The pair-wise related pathway pairs: ");		
		
		Vector<String[]> kwv = comparePathwayNames(kpthv, "K",  wpthv, "W", pw);
		
		Vector<String[]> kcv = comparePathwayNames(kpthv, "K",  cpthv, "C", pw);
		
		Vector<String[]> cwv = comparePathwayNames(cpthv, "C",  wpthv, "W", pw);
		
//		Vector<String[]> kkv = compareInPathwayNames(kpthv, "K", pw);
//		
//		Vector<String[]> wwv = compareInPathwayNames(wpthv, "W", pw);
//		
//		Vector<String[]> ccv = compareInPathwayNames(cpthv, "C", pw);						
		
		Vector<String[]> kkv = comparePathwayNames(kpthv, "K",  kpthv, "K", pw);
		
		Vector<String[]> wwv = comparePathwayNames(wpthv, "W",  wpthv, "W", pw);
		
		Vector<String[]> ccv = comparePathwayNames(cpthv, "C",  cpthv, "C", pw);		
		
        System.out.println("Pairwise Results of pathways that can be merged: \n" +
        		" KEGG-Wiki\tKEGG-BioCyc\twiki-BioCyc\tKEGG-KEGG\tBioCyc-BioCyc \t wiki-wiki \n"
                +kwv.size()+" \t "+kcv.size()+" \t "+cwv.size()+" \t "+kkv.size()+" \t "+ccv.size()+" \t "+wwv.size());
		
        stpw.println("Pairwise Results of pathways that can be merged: \n" +
        		" KEGG-Wiki\tKEGG-BioCyc\twiki-BioCyc\tKEGG-KEGG\tBioCyc-BioCyc \t wiki-wiki \n"
                +kwv.size()+" \t "+kcv.size()+" \t "+cwv.size()+" \t "+kkv.size()+" \t "+ccv.size()+" \t "+wwv.size());
		
        HashMap<String,String> dsjtmp = new HashMap<String,String>();
                
        buildDisjointSet(kwv, dsjtmp);
        
        System.out.println("Pass in kwv (size:"+kwv.size()+")the num of nods disjoint set: "+dsjtmp.size());
        
        buildDisjointSet(kcv, dsjtmp);
        
        System.out.println("Pass in kcv (size:"+kcv.size()+")the num of nods disjoint set: "+dsjtmp.size());
        
        buildDisjointSet(cwv, dsjtmp);
        
        System.out.println("Pass in cwv (size:"+cwv.size()+")the num of nods disjoint set: "+dsjtmp.size());

        buildDisjointSet(kkv, dsjtmp);
        
        System.out.println("Pass in kkv (size:"+kkv.size()+")the num of nods disjoint set: "+dsjtmp.size());
        
        buildDisjointSet(wwv, dsjtmp);
        
        System.out.println("Pass in wwv (size:"+wwv.size()+")the num of nods disjoint set: "+dsjtmp.size());
        
        buildDisjointSet(ccv, dsjtmp);
        
        System.out.println("Pass in ccv (size:"+ccv.size()+") the num of nods disjoint set: "+dsjtmp.size());        

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
    		
    		Vector nodsv = initgroup.get(parent);
    		
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
    		/*
        	String[] rms = pthOriName.split(" ");
        	
        	String un = "";

            for(int x= 0;x<rms.length;x++){
                if(//those conditions are used to modify the integrated pathway names of the pathway groups.  
                   !rms[x].equals("I")&&!rms[x].equals("II")&&!rms[x].equals("V")&&
                   !rms[x].equals("III")&&!rms[x].equals("IV")&&!rms[x].equals("VI")&&
                   !rms[x].contains("(")&&!rms[x].equals(")")
                   ){
                		un = rms[x]+" ";
                	
//                    if(!rms[x].contains("even")&&!rms[x].contains("Type")&&!rms[x].contains("number)")
//                            &&!rms[x].equals("branch)")&&!rms[x].equals("ADP-D-Glucose)")
//                            &&!rms[x].equals("pathway)")){
//                        
//                    }
                    
                }
            }
            */
    		
    		String un = pthOriName;
    		
            if(pathgroup.containsKey(un)){
            	
            	System.out.println("Alert! Something wrong here, two integrated pathway name in different pathway groups.");
            	
            }else{
            	
            	pathgroup.put(un, nodsv);
            	
            }
    		    		    	
    		
    	}
    	
		return pathgroup;
	}


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

    /*
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
            	
            	if (threld >= 0.5) {
                    kwv.add(kws);
                    pw.println(kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
                }
            }

        }
        
        return kwv;
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
                
            	if(threld<0.5) continue;
            	
            	//if (threld >= 0.5) {
                    kwv.add(kws);
                    pw.println(kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
                //}
            }

        }
        
        return kwv;
    }
   
*/
	
//	public static void main(String args[]){
//		
//		String a = "abcdefg";
//		String b = "abc";
//		String c = "efghijk";
//		
//		sequenceAlignment sab = new sequenceAlignment(a, b);
//				
//		sequenceAlignment sac = new sequenceAlignment(a, c);
//		
//	}		

    
	private static Vector<String[]> comparePathwayNames (Vector<String> kv, String kdb,
    		Vector<String> wv, String wdb,PrintWriter pw) throws IOException {
    	
    	Vector<String[]> kwv = new Vector<String[]>();
    	        
        for (int k = 0; k < kv.size(); k++) {        	                     
                    	
            double threld = 0.0;
            
            double aligSc = 0.0;  
            
            double[] param = {0.0, 0.0}; // threld = param[0], aligSc = param[1];
        	
            String pki = kv.get(k);
            
            String pk  = pki.trim();           
        	
    		String[] kws = {"a", "a", "a", "a"};
            
            
            if(kdb.equalsIgnoreCase(wdb)){
            	
                for (int w = k+1; w < wv.size(); w++) {
	            	
                	String wki = wv.get(w);
                	
                	String wk  = wki.trim();
                    
                	
                	seqAlign(pk,kdb,wk,wdb,kws,threld,aligSc);
                	
                	//seqAlign(pk,kdb,wk,wdb,kwv,kws,param);
                }
            	
            	
            }else{
            	
                for (int w = 0; w < wv.size(); w++) {
	            	
                	String wki = wv.get(w);
                	
                	String wk  = wki.trim();
                	
                	seqAlign(pk,kdb,wk,wdb,kws,threld,aligSc);
                        
                }
            	
            }
            
            System.out.println(kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
            
            if (aligSc > kws[0].length() - 1 || aligSc > kws[2].length() - 1 || threld > 0.91) {
                
            	
        		if (threld >= 0.5) {
                    
        			if(!ckmismatch(kws)){
                    	
                    }
        			
                    kwv.add(kws);
                    
                    pw.println(kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
                    
                    //System.out.println(kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
                }
        	}          
        
        }
        
     return kwv;
    
	}





	
//	private static Vector<String[]> buildErrList() {
//		// TODO Auto-generated method stub
//		
//		Vector<String[]> ev = new Vector<String[]>();
//
//		ev.add({ "NOD","Toll"});
//		
//		return ev;
//	}


	private static void seqAlign(String pk, String kdb, String wk, String wdb,
			String[] kws, double threld, double aligSc) {

        //if (wk.length()<2) return;
        
        sequenceAlignment sa = new sequenceAlignment(pk, wk);
        
        if (sa.thres > threld) {//bugs was reported here! Be careful. 
        	aligSc = sa.alsc;
            threld = sa.thres;
            kws[0] = pk;
            kws[1] = kdb;
            kws[2] = wk;
            kws[3] = wdb;            
            System.out.println(kws[0] + "\t" +kws[1] + "\t" + kws[2] + "\t" +kws[3] + "\t" + threld + "\t" + aligSc);
        }   		
	}
		
	//*/
}





