/* Title: Normalizer
 * Function: Java Program to Normalize the information extracted from Markeup Language fies.
 * 
 * Author: Hufeng Zhou
 * Time: July 3rd 2021
 * Version: v2
 */


package sapiens;

import java.io.*;
import java.util.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Set;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class Normalize{	
	
	public static HashMap<String,Vector<String>> mp = new HashMap<String,Vector<String>>(70000);
	
	public static HashMap<String,Vector<String>> testmp = new HashMap<String,Vector<String>>();
		
	public static void buildmap (String orgs)throws IOException{
		
		genecolmp(orgs);
				
		System.out.println("The size of map after retrieve info from HumanCyc genes.col mapping file : " + mp.size());			
		
		
		String uniptf = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"uniprot-homo";
		
		uniptmp(uniptf);
		
		System.out.println("The size of map after retrieve info from Uniprot mapping file : " + mp.size());				
		
		
		//String biomtf = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"mart_export.txt";
		
		//biomtmp(biomtf);
		
		//System.out.println("The size of map after retrieve info from biomart mapping file : " + mp.size());		
		
		//The content in mapping the Entrez to Gene Sybmol needed to write to database.
		
		HashMap<String,Vector<String>> entrzmp = KEGG.ecmp;
		
		System.out.println("The size of map of KO mapping container: "+ entrzmp.size());
		
		mp.putAll(entrzmp);
		
		System.out.println("The size of map after retrieve info from ecmp(KEGG mapping) mapping file : " + mp.size());
		
		
		String HGNCf  = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"NGNC";
						
		hgnc(HGNCf);
		
		System.out.println("The size of map after retrieve info from HGNC mapping file : " + mp.size());
		
		
		String enspro = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"ENStoHGS";
		
		enspmp(enspro);
		
		System.out.println("The size of map after retrieve info from biomart ensembel mapping file : " + mp.size());
		
		String spcf = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"SpecialCorrection";
				
		spmp(spcf);
		
		System.out.println("The size of map after retrieve info from special mapping file : " + mp.size());
		
		HashMap<String,String> widmp = Wiki.wikimp;
		
		prepwikimp(widmp);
		/*
		 * This is the special container contains the main ids of wiki, many of them are fine, few special need to be mapped but the value in this map is not standard ID, conversion is needed. 
		 */								
	}		
	
	private static void enspmp(String enspro) throws IOException {
		
		BufferedReader br = new BufferedReader (new FileReader(enspro));
		
		String ln ;
		
		String[] arr;
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split(",");
			
			if (arr.length<2) continue;
			
			if(arr[0].length()<2||arr[1].length()<2)continue;
			
			String ens = arr[0].trim();
			
			String msb = arr[1].trim();
			
			if(mp.containsKey(ens)){
				
				Vector<String> vv = mp.get(ens);
				
				if(vv.contains(msb)){
					
					vv.add(msb);
				
					mp.put(ens,vv);
				}
				
			}else{
				
				Vector<String> vv = new Vector<String>();
				
				vv.add(msb);
				
				mp.put(ens, vv);
				
			}
			
		}
		
		br.close();
		
	}

	private static void prepwikimp(HashMap<String, String> widmp) {
		
		Set<String> keysets = widmp.keySet();
		
		for(String keys : keysets){
			
			String vals = widmp.get(keys);						
			
			if(!mp.containsKey(keys)&&keys.length()>1){
				
				if(vals.length()<2) continue;
				
				Vector<String> vv = mp.get(vals);
				
				if (vv !=null) {
					
					mp.put(keys,vv);
					
				}
								
			}
			
		}
		
	}

	public HashMap<String,HashMap<String,String>> gemp(HashMap<String,HashMap<String,String>>ige, String db){
		
		mp.remove("ER");
		
		mp.remove("Q6NSA1_HUMAN");
		
		mp.remove("12155");
		
		HashMap<String,HashMap<String,String>> oge = new HashMap<String,HashMap<String,String>>();
		
		//HashMap<String, String> rel = mpRel();
		
		Set pathwayset = ige.keySet();
		
		int oneTomoreCt = 0;
		
		int oneTooneCt= 0;
		
		for (Object pathwayO: pathwayset){
			
			String pathway = pathwayO.toString();
			
			HashMap<String,String> ipgemp = ige.get(pathway);
			
			HashMap<String,String> opgemp = new HashMap<String,String>();
			
			Set pgeset = ipgemp.keySet();
			
			for(Object ipgeO: pgeset){
				
				String ipgek = ipgeO.toString();
				
				String opgek = "";
				
				if(mp.containsKey(ipgek)){
					
					Vector<String> vk = mp.get(ipgek);
					
					if(vk.size()>1){
						oneTomoreCt++;
					}else{
						oneTooneCt++;
					}
					
					for(int i = 0; i<vk.size(); i++){
						
//						if(vk.get(i).contains("/")){
//							
//							String[] tmp = vk.get(i).split("\\/");
//							
//							String opgek0 = tmp[1];
//							
//						}
						
						String opgekPre = vk.get(i).trim();
						
						opgek = opgekPre.replaceAll("\\?", "").replaceAll(" ", "").replaceAll("`","").replaceAll("\\+","");//.toUpperCase()
						
						
						
						if(!opgemp.containsKey(opgek)&&opgek.length()>1){
							
							opgemp.put(opgek,db);
							
						}
						
					}
					
				}else{
					
					opgek = ipgek;
					
					if(opgek.length()>1) {
					
						opgemp.put(opgek,db);
						
					}					
															
				}
				
			}
			
			oge.put(pathway, opgemp);
			
		}
		
		System.out.println( "In mapping genes in "+db+"\n"+
				"The number of genes have one to more mapping: "+oneTomoreCt+
				"\nThe number of genes have one to one mapping:"+oneTooneCt
				);
		
		return oge;
		
	}
	
	public HashMap<String,HashMap<String,String>> gpmp(HashMap<String,HashMap<String,String>>igp)throws Exception{
		
		mp.remove("ER");
		
		mp.remove("Q6NSA1_HUMAN");
		
		mp.remove("12155");
		
		HashMap<String,HashMap<String,String>> gp = new HashMap<String,HashMap<String,String>>();
				
		HashMap<String, String> rel = mpRel();
		
		Set<String> ipgpset = igp.keySet();
		
		int oneTomoreCt = 0;
		
		int oneTooneCt= 0;
		
		for(String pwy : ipgpset){
			
			//String pwy = pwyO.toString();
			
			HashMap<String,String> igpmp  = new HashMap<String,String>();
			
			HashMap<String,String> opgpmp = new HashMap<String,String>();
			
			igpmp = igp.get(pwy);
			
			Set<String> igpset = igpmp.keySet();						
			
			for(String igps : igpset){				
				
				String irels = igpmp.get(igps); 
				
				String[] ge = igps.split("\t");
				
				Vector<String> va = new Vector<String>();
				
				Vector<String> vb = new Vector<String>();
				
				if(mp.containsKey(ge[0])&& mp.containsKey(ge[1])){
					
					va = mp.get(ge[0]);
					
					vb = mp.get(ge[1]);										
					
				}else if(mp.containsKey(ge[0])&& !mp.containsKey(ge[1])){
					
					va = mp.get(ge[0]);
					
					vb.add(ge[1]);
					
					
				}else if(!mp.containsKey(ge[0])&& mp.containsKey(ge[1])){
					
					va.add(ge[0]);
					
					vb = mp.get(ge[1]);	
					
					
					
				}else if(!mp.containsKey(ge[0])&&!mp.containsKey(ge[1])){					
					
					va.add(ge[0]);
					
					vb.add(ge[1]);
					
				}
				
				if(va.size()>1){
					
					oneTomoreCt++;
					
					if(! testmp.containsKey(ge[0])){
						
						testmp.put(ge[0], va);
						
					}
					
				}else{
					
					oneTooneCt++;
					
				}
				
				if(vb.size()>1){
					
					oneTomoreCt++;
					
					if(! testmp.containsKey(ge[1])){
						
						testmp.put(ge[1], va);
						
					}
					
				}else{
					
					oneTooneCt++;
					
				}
				
				
				String orels = rel.get(irels);							
				
				for(int i = 0; i< va.size(); i++){
					
					for(int j = 0; j< vb.size();j++){
						
						String oga = "";
						
						String ogb = "";
						
						if(va.get(i).compareTo(vb.get(j))>=0){
							
							oga = va.get(i);
							
							ogb = vb.get(j);
							
						}else{
							
							oga = vb.get(j);
							
							ogb = va.get(i);
														
						}
						
						String oa = oga.replaceAll("\\?", "").replaceAll(" ", "").replaceAll("`","").replaceAll("\\+","");//.toUpperCase()
						
						String ob = ogb.replaceAll("\\?", "").replaceAll(" ", "").replaceAll("`","").replaceAll("\\+","");//.toUpperCase()
						
						if(oa.length()>1&&ob.length()>1){
							
							StoreGenePairs(oa.trim(),ob.trim(), orels,pwy,gp);
							
						}
						
					}
				}
			}
			
			
		}
		
		System.out.println(
				"The number of genes in gene-pairs have one to more mapping: "+oneTomoreCt+
				"\nThe number of genes in gene-pairs have one to one mapping:"+oneTooneCt
				);
		
		return gp;
		
	}	
		
    private static void StoreGenePairs(String geneA, String geneB, String rel,String pathway, 
    		HashMap<String,HashMap<String,String>> phgp) {
		
    	HashMap<String, String> gp = new HashMap<String, String>();
    	
    	String ga, gb;
    	
    	if(geneA.compareTo(geneB)>=0){    		
    		
    		ga = geneA;
    		
    		gb = geneB;
    		
    	}else{
    		
    		ga = geneB;
    		
    		gb = geneA;
    		
    	}
    	
    	if(phgp.containsKey(pathway)){
    		
    		gp = phgp.get(pathway);
    		
    		if(!gp.containsKey(ga+"\t"+gb)){
    			
    			gp.put(ga+"\t"+gb, rel);
    			
    			phgp.put(pathway, gp);
    			
    		}    		    	
    		
    	}else{
    		
    		gp.put(ga+"\t"+gb, rel);
    		
    		phgp.put(pathway, gp);
    		
    	}
		
	}  	
	
	public static HashMap<String,String> mpRel(){

        HashMap<String,String> relmp = new HashMap<String,String>();
        
        relmp.put("SEQUENTIAL_CATALYSIS", "ECrel");
        relmp.put("ECrel","ECrel");
        relmp.put("IN_SAME_COMPONENT", "GPrel");//PPrel in V1.0
        //relmp.put("COMPONENT_OF", "GPrel");//no, PPrel in V1.0 
        relmp.put("CO_CONTROL", "PPrel");
        relmp.put("INTERACTS_WITH", "PPrel");//no
        relmp.put("PPrel","PPrel");
        relmp.put("PCrel","PPrel");// Correct Error.
        relmp.put("GErel","GErel");//only in KEGG
        relmp.put("maplink","maplink");//no
        relmp.put("groupRel","GPrel");//PPrel in V1.0
        relmp.put("graphRel","PPrel");
        relmp.put("GPrel","GPrel"); //newly added.
        
        return relmp;
		
	}

	
	public static void genecolmp(String orgs) throws IOException{
		
		String gecol   = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"genes.col";
		
		String[] arr;
		
		String ln;
				
		BufferedReader br = new BufferedReader(new FileReader(gecol));
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split("\t");
			
			if(arr.length<3) continue;
			
			if(arr[1].length()<2) continue;		
			
			Vector<String> kv = new Vector<String>();			
			
			//int[] indx = {0,3,7,8,9,10};
	
			if(arr[0].length()>1){
				
				kv.add(arr[0].trim());
				
			}
			
//			if(arr[3].length()>1){
//				
//				kv.add(arr[3]);
//				
//			}			
			
			
//			for(int i = 7;i<arr.length && i<11;i++){
//								
//				if(arr[i].length()>1){
//					
//					kv.add(arr[i]);
//					
//				}
//				
//			}
			
//			if(arr[2].length()>1 && arr[2].length()<20){
//				
//				kv.add(arr[2]);
//				
//			}
			
//			if(arr[1].length()>1){
//				
//				Vector<String> tmp = new Vector<String>();
//				
//				tmp.add(arr[1]);
//				
//				mp.put(arr[1],tmp);
//				
//			}	
			
			continueBuildmp(kv, arr[1].trim());
			
//			for(int i = 0; i<kv.size();i++){
//				
//				Vector<String> vv = new Vector<String>();
//				
//				if(mp.containsKey(kv.get(i))){
//					
//					vv = mp.get(kv.get(i));
//					
//					if(!vv.contains(arr[1])){
//						
//						vv.add(arr[1]);
//						
//						mp.put(kv.get(i),vv);
//					}
//					
//				}else{
//					
//					vv.add(arr[1]);
//					
//					mp.put(kv.get(i),vv);
//					
//				}
//				
//			}
			
		}		
		
		br.close();
		
	}
		
	private static void uniptmp(String uniptf) throws IOException {
        
		BufferedReader mpbr = new BufferedReader(new FileReader(uniptf));
		
		mpbr.readLine();
		
		String lnm;

		String[] arr;

		while((lnm=mpbr.readLine())!=null){		    	
			
	        arr = lnm.split("\t");			     
	        
	        if(arr[4].length()<2) continue;	        

	        Vector<String> oidv = new Vector<String>();
				
			String[] ge = arr[4].split(" ");
			
			for(int i = 1;i< ge.length;i++) {

			    if (ge[i].length()>0) oidv.add(ge[i].trim());

			}		
						
			
			if(arr[0].length()>1){
			
				oidv.add(arr[0].trim());
				
			}
		        
			if(arr[1].length()>1){
			
				oidv.add(arr[1].trim());
			
			}
		        

			//Pattern ep = Pattern.compile("EC (\\d{1,4}\\.\\d{0,4}\\.\\d{0,4}\\.[\\d|-]{0,4})");
			
			Pattern ep = Pattern.compile("EC (\\d{1,4}\\.\\d{0,4}\\.\\d{0,4}\\.[\\d|-]{0,4})");

			Matcher em = ep.matcher(arr[3]);

			while(em.find()){							

				if(em.group(1).length()<1) break;

				Vector<String> vec = new Vector<String>();

				oidv.add(em.group().trim());
				
				oidv.add(em.group(1).trim());
				
				//System.out.println(em.group());
				
				//System.out.println(em.group(1));
				
			}			

	        String kid = ge[0].trim(); //kid is the unified gene id
	        
	        if(kid.length()<1) continue;

			//Vector<String> kidv = new Vector<String>();
			
	        continueBuildmp(oidv, kid.trim());
	        
//	        for(int j = 0; j< oidv.size();j++){
//
//	                String oid = oidv.get(j); // oid is the other gene names.
//
//	                //gene_mapping	
//	                
//	                Vector<String> kidv = new Vector<String>();
//
//	                if(mp.containsKey(oid)){	                		                	
//	                	
//						kidv = mp.get(oid);                  
//						
//						if(!kidv.contains(kid)){
//		
//							kidv.add(kid);
//							
//							mp.put(oid,kidv);
//							
//						} 										
//	                        
//	                }else{	
//	                				
//	                        kidv.add(kid);
//	                        
//	                        mp.put(oid,kidv);
//	                        
//	                }            
//				
//	        }

		}

		mpbr.close();
    }	
	
	private static void biomtmp(String biomtf) throws IOException { //, HashMap<String,Vector<String>>mp
		
		String[] arr;
		
		String ln;		
		
		BufferedReader br = new BufferedReader(new FileReader(biomtf));
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split(",");
			
			if(arr.length<4) continue;
			
			if(arr[1].length()<2) continue;
			
			Vector<String> kv = new Vector<String>();				
			
			if(arr[0].length()>1){
				
				kv.add(arr[0].trim());
				
			}
			
			if(arr[1].length()>1){
				
				Vector<String> tmp = new Vector<String>();
				
				tmp.add(arr[1]);
				
				mp.put(arr[1].trim(),tmp);
				
			}
			
			if(arr[2].length()>1){
				
				kv.add(arr[2].trim());
				
			}
			
			if(arr[3].length()>1){
				
				kv.add(arr[3].trim());
				
			}
			
			continueBuildmp(kv, arr[1].trim());			
			
		}
		
		
		br.close();				
	}	
	
	private static void continueBuildmp(Vector<String> kv, String sid) {
		
		if(sid.contains("'")||sid.contains("ase")) return;
		
		for(int i = 0; i<kv.size();i++){
			
			if(kv.get(i).length()<2) return;
			
			Vector<String> vv = new Vector<String>();
			
			if(mp.containsKey(kv.get(i))){
				
				vv = mp.get(kv.get(i));
				
				if(!vv.contains(sid)){
					
					vv.add(sid);
					
					mp.put(kv.get(i),vv);
				}
				
			}else{
				
				vv.add(sid);
				
				mp.put(kv.get(i),vv);
				
			}
			
		}		
		
	}

	private static void hgnc(String ncbif) throws IOException {
		
		String[] arr;
		
		String ln;
				
		BufferedReader br = new BufferedReader(new FileReader(ncbif));
		
		br.readLine();
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split("\t");
			
			if(arr.length<8) continue; 
			
			if(arr[1].length()<2) continue;
			
			Vector<String> kv = new Vector<String>();									
			
			if(arr[0].length()>1){
				
				//String[] ar0 = arr[0].split("\\:");
				
				kv.add(arr[0].trim());
				
				//System.out.println(ar0[1]);
				
			}
			
			if(arr[6].length()>1){
				
				if(arr[6].contains(",")){
					
					String[] tmp = arr[6].split(",");
					
					for(int s = 0; s<tmp.length;s++){
						
						kv.add(tmp[s].trim());
						
					}					
					
				}else{
					
					kv.add(arr[6].trim());
					
				}
				
			}	
			
			if(arr[7].length()>1){
				
				if(arr[7].contains(",")){
					
					String[] tmp = arr[7].split(",");
					
					for(int s = 0; s<tmp.length;s++){
						
						kv.add(tmp[s].trim());
						
					}					
					
				}else{
					
					kv.add(arr[7].trim());
					
				}
				
			}			
		
			/*
			 * This part previous symbol introduce many noise from the HGNC symbol mapping file.
			 */
//			if(arr[3].length()>1){
//				
//				if(arr[3].contains(",")){
//					
//					String[] tmp = arr[3].split(",");
//					
//					for(int s = 0; s<tmp.length;s++){
//						
//						kv.add(tmp[s].trim());
//						
//					}					
//					
//				}else{
//					
//					kv.add(arr[3].trim());
//					
//				}				
//										
//			}
			/*
			 * This is used to get rid of errors from synomonys
			 */
			
			
//			if(arr[4].length()>1){
//				
//				if(arr[4].contains(",")){
//					
//					String[] tmp = arr[4].split(",");
//					
//					for(int s = 0; s<tmp.length;s++){
//						
//						kv.add(tmp[s].trim());
//						
//					}					
//					
//				}else{
//					
//					kv.add(arr[4].trim());
//					
//				}									
//				
//			}				
			
			if(arr[1].length()>1){
				
				Vector<String> tmp = new Vector<String>();
				
				tmp.add(arr[1].trim());
				
				mp.put(arr[1].trim(),tmp);
				
			}
			
			continueBuildmp(kv, arr[1].trim());
			
//			for(int i = 0; i<kv.size();i++){
//				
//				Vector<String> vv = new Vector<String>();
//				
//				if(mp.containsKey(kv.get(i))){
//					
//					vv = mp.get(kv.get(i));
//					
//					if(!vv.contains(arr[1])){
//						
//						vv.add(arr[1]);
//						
//						mp.put(kv.get(i),vv);
//					}
//					
//				}else{
//					
//					vv.add(arr[1]);
//					
//					mp.put(kv.get(i),vv);
//					
//				}
//				
//			}
			
		}		
		
		br.close();		
		
	}	

	//This is just a quick and dirty way to handle some pathway IDs that failed to map.
	private static void spmp(String spcf) throws IOException{
		
		BufferedReader br = new BufferedReader(new FileReader(spcf));
		
		String ln;
		
		String[] arr;
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split("=");
			
			spf(arr[0],arr[1]);
			
		}
		
		br.close();
		
	}

	private static void mapUpcase() {
		
		String[] uppercase = {
//				"","","","","",
//				"","","","","",
//				"","","","","",
//				"","","","","",
//				"","","","","",
				};
		
		for(int i = 0; i<uppercase.length;i++){
		
			String inputID = uppercase[i];
			
			String outputID = inputID.substring(0,1)+inputID.substring(1,inputID.length()).toLowerCase();
			
			spf(inputID,outputID);
			
		}
		
	}

	private static void spf(String iName, String oName) { //iName is input , oName is output standard name.
		
		Vector<String> v = new Vector<String>();
		
		v.add(iName);
		
		continueBuildmp(v,oName);
		
		
		
	}	
	

	public HashMap<String, HashMap<String, Vector<String>>> groupmp(
			HashMap<String, HashMap<String, Vector<String>>> wprgro) {

		
		HashMap<String,HashMap<String,Vector<String>>> gro = new HashMap<String,HashMap<String,Vector<String>>>();
		
		//HashMap<String, String> rel = mpRel();
		
		Set<String> pathwayset = wprgro.keySet();
				
		for (String pathway: pathwayset){			
			
			HashMap<String,Vector<String>> igro = wprgro.get(pathway);
			
			HashMap<String,Vector<String>> ogro = new HashMap<String,Vector<String>>();
			
			Set<String> groidset = igro.keySet();
			
			for(String groids: groidset){
				
				Vector<String> inv = igro.get(groids);
				
				Vector<String> otv = new Vector<String>();
				
				for(int x = 0; x < inv.size();x++){
					
					String inge = inv.get(x);
					
					String otge = "";
					
					if(mp.containsKey(inge)){
						
						Vector<String> vk = mp.get(inge);
						
						
						for(int i = 0; i<vk.size(); i++){
							
							
							String opgekPre = vk.get(i);
							
							otge = opgekPre.replaceAll("\\?", "").replaceAll(" ", "").replaceAll("`","").replaceAll("\\+","");//.toUpperCase()
							
							
							
							if(!otv.contains(otge)&&otge.length()>1){
								
								otv.add(otge);
								
							}
							
						}
						
					}else{
						
						otge = inge;
						
						if(!otv.contains(otge)&&otge.length()>1){
							
							otv.add(otge);
							
						}				
																
					}
					
				}
				
				//ogro.put("group", value);
				
            	if(!ogro.containsKey("Group"+groids)){
            		
            		ogro.put("Group"+groids,otv);
            		
            	}else{
            		
            		System.out.println("Alert! It seems some group ids are the same!");
            		
            		otv.addAll(ogro.get("Group"+groids));
            		
            		ogro.put("Group"+groids,otv);
            		
            	}
					
					
			}
				
			gro.put(pathway, ogro);
			
		}
		
		System.out.println( "Number of groups: "+gro.size());
		
		return gro;
	}		
	
	public static void main(String args[])throws IOException{
		
		buildmap("sapiens");
		
	}
}


