package tuberculosis;

import java.io.*;
import java.util.*;
import java.util.regex.*;


public class Normalize{	
	
	public static HashMap<String,Vector<String>> mp = new HashMap<String,Vector<String>>(30000);
	
	public static HashMap<String,Vector<String>> testmp = new HashMap<String,Vector<String>>();
	
	public static void main(String args[]) throws IOException{
		
		System.out.println("Now using the normaliza main class to debug...");
		
		buildmap("tuberculosis");
		
	}
		
	public static void buildmap (String orgs)throws IOException{
		
    	
		String ecmpf = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"mtu_enzyme.list";
    	
		EcSybMap(ecmpf);
		
		HashMap<String,Vector<String>> entrzmp = KEGG.ecmp;
		
		System.out.println("The size of map of KO mapping container: "+ entrzmp.size());
		
		mp.putAll(entrzmp);
		
		System.out.println("The size of map after retrieve info from ecmp(KEGG mapping) mapping file : " + mp.size());		
		
		String simplempf = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"mtuRvGN";
		
		simplemp(simplempf);
		
		System.out.println("The size of map after retrieve info from simplemp mapping file : " + mp.size());		
			
//		String geneNameMp = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"SGD_features.tab";
//		
//		MapGeNtoSysN(geneNameMp);
//		
//		System.out.println("The size of map after retrieve info from SGD Symbol name to systematic name mapping file : " + mp.size());						
//		
//		
//    	String ncbi 	= orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"NCBImapping";    	    	
//   	
//    	EntrSybNCBI(ncbi);		
//    	
//    	System.out.println("Size of mapping file after ncbi mapping file: "+mp.size());
//    	
//    	
//    	String gecol   = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"genes.col";
//    	
//    	genecolmp(gecol);
//    	
//    	System.out.println("The size of map after retrieve info from BioCyc gene mapping file : " + mp.size());
//    	
//		
//		String uniportf = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"sce_uniprot.list";
//		
//		uniport(uniportf);
//		
//		System.out.println("The size of map after retrieve info from simplemp mapping file : " + mp.size());
//    	
    	
		String spcf = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"SpecialCorrection";
				
		spmp(spcf);
		
		System.out.println("The size of map after retrieve info from special mapping file : " + mp.size());
		
		HashMap<String,String> widmp = Wiki.wikimp;
		
		System.out.println("The size of wikimap is: "+widmp.size());
		
		prepwikimp(widmp);
		
		/*
		 * This is the special container contains the main ids of wiki, many of them are fine, few special need to be mapped but the value in this map is not standard ID, conversion is needed. 
		 */								
	}		
	
	private static void EcSybMap(String kompf) throws IOException{
		// this method is used to include more information to the ecmp container, koid as key, GeneSymbol as value in Vector.
		
		String ln ;
		
		String[] arr;
		
		BufferedReader br = new BufferedReader(new FileReader(kompf));
		
		while((ln= br.readLine())!=null){
			
			arr = ln.split("\t");
			 
			if(arr.length<2||arr[0].length()<2||arr[1].length()<2){
				continue;
			}
			
			//System.out.println(ln);
			
			Vector<String> tv  = new Vector<String>();
			
			if(mp.containsKey(arr[1].trim())){
				
				tv = mp.get(arr[1].trim());
				
				if(!tv.contains(arr[0].trim())){
					
					tv.add(arr[0].trim());
					
					mp.put(arr[1].trim(), tv);
					
				}								
				
			}else{
				
				Vector<String> vv = new Vector<String>();
				
				vv.add(arr[0].trim());
				
				mp.put(arr[1].trim(), vv);
				
			}		
			
			if(arr[0].length()>1){
				
				Vector<String> tmp = new Vector<String>();
				
				tmp.add(arr[0].trim());
				
				mp.put(arr[0].trim(),tmp);
				
			}
		}
		
		
		br.close();
		
	}
	
	private static void uniport(String uniportf) throws IOException{
		String[] arr;
		
		String ln;		
		
		BufferedReader br = new BufferedReader(new FileReader(uniportf));
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split("\t");
			
			if(arr.length<2) continue;
			
			if(arr[1].length()<2) continue;
			
			Vector<String> kv = new Vector<String>();				
			
			if(arr[1].length()>1){
				
				kv.add(arr[1].trim());
				
			}			
			
			
			
			if(kv.size()>0){
				
				//System.out.println(arr[0]+"\t"+arr[1]);
				
				continueBuildmp(kv, arr[0].trim());	
				
			}
			
			
			
			
		}
		
		
		br.close();		
		
	}

	public static void genecolmp(String gecol) throws IOException{		
		
		String[] arr;
		
		String ln;
				
		BufferedReader br = new BufferedReader(new FileReader(gecol));
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split("\t");
			
			if(arr.length<8) continue;
			
			if(arr[0].length() <2 ) continue;		
			
			Vector<String> kv = new Vector<String>();			
			
	
			if(arr[0].contains("G") || arr[0].contains("Q")){
				
				kv.add(arr[0].trim());
								
				
			}
			
			if(arr[7].contains("Y")&&arr[7].length()<15&&kv.size()>0){
				
				continueBuildmp(kv, arr[7].trim());	
				
				//System.out.println(arr[0].trim()+"\t"+arr[7]);
				
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
		
	private static void MapGeNtoSysN(String input) throws IOException {
        
		BufferedReader mpbr = new BufferedReader(new FileReader(input));
		
		mpbr.readLine();
		
		String lnm;

		String[] arr;

		while((lnm=mpbr.readLine())!=null){		    	
			
	        arr = lnm.split("\t");			     
	        
	        if(arr.length<5) continue; 
	        
	        //if(arr[4].length()<2) continue;	        

	        Vector<String> oidv = new Vector<String>();
						
						
			
			if(arr[4].length()>1){
			
				oidv.add(arr[4].trim());
				
			}
		        
			if(arr[4].length()>3 && arr[3].contains("Y")){
				
		        String kid = arr[3].trim(); //kid is the unified gene id
		        
		        //System.out.println(kid+"\t"+arr[4]);
		        
		        continueBuildmp(oidv, kid.trim());
			}


		}

		mpbr.close();
    }	
	
	

	private static void EntrSybNCBI(String ncbi) throws IOException {
		// TODO Auto-generated method stub
		String[] arr;
		
		String ln;		
		
		BufferedReader br = new BufferedReader(new FileReader(ncbi));
		
		String symbo, entrz;
		
		symbo= "";
		
		entrz = "";
		
		int indx = 0;
		
		while(br.ready()){//(ln = br.readLine())!=null
			
			ln = br.readLine();
			
			indx++;
			
			if(ln.length()<1){
				
				if(symbo==""||entrz =="") continue;
				
				Vector<String> kv = new Vector<String>();
				
				kv.add(entrz.trim());
				
				continueBuildmp(kv, symbo.trim());
				
				//System.out.println(symbo+"\t"+entrz);
				
				symbo= "";
				
				entrz = "";
				
				
			}
			
			if(ln.contains("Other Aliases")){
				
				//System.out.println(ln);
				
				Pattern sbp = Pattern.compile("Other Aliases\\: (Y[\\w|-]{5,10})[,|\\n]{1}");
				
				Matcher sbm = sbp.matcher(ln);
				
				while(sbm.find()){
					
					symbo = sbm.group(1);
					
					//System.out.println(symbo);
					
				}
				
			}
			
			if(ln.contains("ID:")){
				
				Pattern idp = Pattern.compile("ID\\:(.*)");
				
				Matcher idm = idp.matcher(ln);
				
				while(idm.find()){
					
					entrz = idm.group(1);
					
				}
				
			}
			
			
		}
		
		///*
		System.out.println(//"The ncbi mapping file have just finish construction, size: "+mp.size()+
				"\nAnd how many lines of the mapping file: "+ indx);
		//*/
	}
	

	
	private static void simplemp(String simplempf) throws IOException { //, HashMap<String,Vector<String>>mp

	    String ln;
	    String[] arr;

	    BufferedReader br = new BufferedReader(new FileReader(simplempf));
	    
	    while((ln = br.readLine())!=null){
	        
	    	arr = ln.split("\t");
	        
	        if(arr[1].length()>1){
	            
	        	Vector<String> grv = new Vector<String>();
	            
	            if(!arr[1].contains(",")){
	                
	            	if(!mp.containsKey(arr[1].trim())){
	                    
	                	grv.add(arr[0].trim());
	                    mp.put(arr[1].trim(),grv);
	                    
	                }else if(mp.containsKey(arr[0].trim())){
	                    
	                	if(!grv.contains(arr[1].trim())){
	                        grv.add(arr[0].trim());
	                        mp.put(arr[1].trim(), grv);
	                    }
	                    
	                    //System.out.println("Bugs:More than one gene name exist: Key>"+arr[1]+"  Value>"+arr[0]);
	                }
	            }else if(arr[1].contains(",")){
	                String[] gk = arr[1].split(",");
	                    for(int k = 0; k<gk.length;k++){
	                    if(!mp.containsKey(gk[k].trim())){
	                        grv.add(arr[0].trim());
	                        mp.put(gk[k].trim(),grv);
	                    }else if(mp.containsKey(gk[k].trim())){
	                        if(!grv.contains(arr[0].trim())){
	                            grv.add(arr[0].trim());
	                            mp.put(gk[k].trim(), grv);
	                        }
	                        
	                        //System.out.println("Bugs:More than one gene name exist: Key>"+arr[1]+"  Value>"+arr[0]);
	                    }
	                }
	            }

	        }
	        
			if(arr[0].length()>1){
				
				Vector<String> tmp = new Vector<String>();
				
				tmp.add(arr[0].trim());
				
				mp.put(arr[0].trim(),tmp);
				
			}
	    }
	    br.close();
		
		
//		String[] arr;
//		
//		String ln;		
//		
//		BufferedReader br = new BufferedReader(new FileReader(simplempf));
//		
//		while((ln = br.readLine())!=null){
//			
//			arr = ln.split("\t");
//			
//			if(arr.length<2) continue;
//			
//			if(arr[1].length()<2) continue;
//			
//			if(arr[1].length()>1){
//				
//				Vector<String> tmp = new Vector<String>();
//				
//				tmp.add(arr[1]);
//				
//				mp.put(arr[1].trim(),tmp);
//				
//			}
//			
//			Vector<String> kv = new Vector<String>();				
//			
//			if(arr[0].length()>1){
//				
//				kv.add(arr[0].trim());
//				
//				continueBuildmp(kv, arr[1].trim());
//				
//				//System.out.println(arr[1]+"\t"+arr[0]);
//				
//			}			
//			
//						
//			
//		}
//		
//		
//		br.close();				
	}	
	
	private static void continueBuildmp(Vector<String> kv, String sid) {
		
		//if(sid.contains("'")||sid.contains("ase")) return;
		
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
	
}



