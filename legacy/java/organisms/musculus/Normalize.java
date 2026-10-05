package musculus;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class Normalize{	
	
	public static HashMap<String,Vector<String>> mp = new HashMap<String,Vector<String>>();
	
	public static HashMap<String,Vector<String>> testmp = new HashMap<String,Vector<String>>();
	
	public static void buildmap (String orgs)throws IOException{
		
		String uniptf = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"uniprot-Mus";
		
		String biomtf = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"BioMart-MousefMGI";
		
		String vega   = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"MGIVEGA";
		
		String ncbif  = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"GeneBankID";//
		
		String enspro = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"EnsProteinMus";//ENStoHGS
		
		enspmp(enspro);
		
		uniptmp(uniptf);
		
		System.out.println("The size of map after retrieve info from Uniprot mapping file : " + mp.size());
		
		biomtmp(biomtf);
		
		System.out.println("The size of map after retrieve info from BioMart mapping file : " + mp.size());
		
		vegamp(vega);
		
		System.out.println("The size of map after retrieve info from MGI-VEGA mapping file : " + mp.size());
		
		ncbi(ncbif);
		
		System.out.println("The size of map after retrieve info from NCBI mapping file : " + mp.size());		
		
		String spcf = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"SpecialCorrection";
		
		
		//The content in mapping the Entrez to Gene Sybmol needed to write to database.
		
		HashMap<String,Vector<String>> entrzmp = KEGG.ecmp;
		
		System.out.println("The size of map of KO mapping container: "+ entrzmp.size());
				
		mp.putAll(entrzmp);		
		
		System.out.println("The size of map after retrieve info from KO mapping file : " + mp.size());
		
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
						
						String opgekPre = vk.get(i);
						
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
				
				if(!rel.containsKey(irels))System.out.println("There are some relationships can't be matched"+irels);
				
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
						
						String oa = oga.replaceAll("\\?", "").replaceAll(" ", "").replaceAll("`","").replaceAll("\\+","").trim();//.toUpperCase()
						
						String ob = ogb.replaceAll("\\?", "").replaceAll(" ", "").replaceAll("`","").replaceAll("\\+","").trim();//.toUpperCase()
						
						if(oa.length()>1&&ob.length()>1){
							
							StoreGenePairs(oa,ob, orels,pwy,gp);
							
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

			    if (ge[i].length()>0) oidv.add(ge[i]);

			}		
						
			
			if(arr[0].length()>1){
			
				oidv.add(arr[0].trim());
				
			}
		        
			if(arr[1].length()>1){
			
				oidv.add(arr[1].trim());
			
			}
		        

			Pattern ep = Pattern.compile("EC (\\d{1,4}\\.\\d{0,4}\\.\\d{0,4}\\.[\\d|-]{0,4})");

			Matcher em = ep.matcher(arr[3]);

			while(em.find()){
			//System.out.println(em.group(1));

				if(em.group(1).length()<1) break;

				Vector<String> vec = new Vector<String>();

				oidv.add(em.group());
								
				oidv.add(em.group(1));				
				
			}			

	        String kid = ge[0]; //kid is the unified gene id
	        
	        if(kid.length()<1) continue;
	        
	        continueBuildmp(oidv, kid.trim());

		}

		mpbr.close();
    }	
	
	private static void biomtmp(String biomtf) throws IOException { //, HashMap<String,Vector<String>>mp
		
		String[] arr;
		
		String ln;		
		
		BufferedReader br = new BufferedReader(new FileReader(biomtf));
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split(",");
			
			if(arr[1].length()<2) continue;
			
			Vector<String> kv = new Vector<String>();
			
			Vector<String> vv = new Vector<String>();
			
			if(arr[0].length()>1){
				
				kv.add(arr[0].trim());
				
			}
			
			if(arr[1].length()>1){
				
				Vector<String> tmp = new Vector<String>();
				
				tmp.add(arr[1].trim());
				
				mp.put(arr[1].trim(),tmp);
				
			}
			
//			if(arr[2].length()>1){
//				
//				kv.add(arr[2].trim());
//				
//			}
			
			if(arr[3].length()>1){
				
				kv.add(arr[3].trim());
				
			}
			
			continueBuildmp(kv, arr[1].trim());
			
		}
		
		
		br.close();				
	}	
	
	private static void vegamp(String vega) throws IOException{
	
		
		String[] arr;
		
		String ln;
				
		BufferedReader br = new BufferedReader(new FileReader(vega));
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split("\t");
			
			if(arr[1].length()<2) continue;
			
			Vector<String> kv = new Vector<String>();
			
			//Vector<String> vv = new Vector<String>();
			
			if(arr[5].length()>1){
				
				kv.add(arr[0].trim());
				
			}
			
			if(arr[1].length()>1){
				
				Vector<String> tmp = new Vector<String>();
				
				tmp.add(arr[1].trim());
				
				mp.put(arr[1].trim(),tmp);
				
			}						
			
			continueBuildmp(kv, arr[1].trim());						
			
		}		
		
		br.close();	
		
		
	}
	
	private static void ncbi(String ncbif) throws IOException {
		
		String[] arr;
		
		String ln;
				
		BufferedReader br = new BufferedReader(new FileReader(ncbif));
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split(",");
			
			if(arr.length<3) continue; 
			
			if(arr[2].length()<2) continue;
			
			Vector<String> kv = new Vector<String>();
			
			//Vector<String> vv = new Vector<String>();
			
			
			if(arr[0].length()>1){
				
				kv.add(arr[0].trim());
				
			}
			
			
			if(arr[1].length()>1){
				
				kv.add(arr[1].trim());
				
			}
			
			continueBuildmp(kv, arr[2].trim());								
			
		}		
		
		br.close();		
		
	}	
	
	//This is just a quick and dirty way to handle some pathway IDs that failed to map.
	
	
	private static void mapUpcase() {
		
		String[] uppercase = {"CYP2A13", "MTHFD1L", "MTHFD1", "MTHFS", "MTHFD2", "TYMS", "MTHFR", "ATIC", 
				"MAT1A", "TCN II", "SHMT1", "SHMT2", "FTCD", "GART", "MTRR", "FOLH1", "DHFR", "MTR", "CYP2B6",
				"BHMT", "DNMT1", "MAT2B", "MTFMT", "CYP3A4", "G1P3", "OAS1", "SAP18", "CCNB1", "CYP2D6", 
				"ADH1C", "ADH1B", "ADH6", "ADH1A","RPS6KB1","SKP1A","KPNA2","UBCH5B","SUMO1","COXIII",
				"COXII","CYTB","NDUFB1","ND6","CYP4A11","UGT2B7","HIST2H4","CYP4A11","TNFRSF10B","Cdc14A","Cdc14B",
				"CYP2J2","CYP2C9","CYP2C8","CYP3A5","CYP3A7","CYP2J2","CYP2C19","CYP2C18","CYP2G1P","RBPSUHL",
				"AHCY","AKR1B10","ARSC","ARSD","ARSE","GAL3ST2","IGHA2","GRM6","GSTM1","KCNQ1","MAP3K7","MMP26","V1RC9",
				"GJA10","GLYATL1","GLYATL2","GNG5","GST3","GSTA5","GSTA4","GSTA3","GSTA2","GSTA1","IGHA1",
				"GSTM5","GSTM4","GSTP1","HADHB","HIST2H3C","HIST3H3","HS3ST3A1","HS3ST3B1","HSPA1A","ICAM3",
				"IRAK4","IRAK3","IL3","IKBKAP","IGHM","UGT2B4","KLF15","MAGEA1","MAP2K1","MAP3K1","MAP3K14",
				"MAPK1","MAPK3","MAPK8","MAPT","MAPT","NCL","NFKBIB","NFKBIA","NFKB1","neuroD4","neuroD1","NCOA4",
				"NRF2B1","PCNA","PDK1","PELI1","PGM5","PI3K","PRKACG","PPP2R3B","PLCG1","PLAGL1","PILRB","PIK3R1",
				"PRKCG","PRKCZ","PSENEN","PTPN11","PTPN3","RELA","SelH","SCTR","RUNX1","RPL6","RPL30","RPA4","UGT2B28",
				"SelK","SelM","SelO","SelS","SelV","SelT","UGT1A9","SepW1","SepX1","SOCS2","SQSTM1","SRY","STAT3",
				"SULT1C2","SULT2A1","TOLLIP","UGT1A10","UGT1A7","SULT1C3","UGT2A2","UGT2B11","UGT2B15","UGT2B17",
				"SULT1C4","CYP4F2","CYP4F11","CYP4F12","CYP4A22","CYP4Z1","ENAH","EPPK1","NCF1","NAT1","NAT2",
				"NAT8","NCF1","PI5K","SEPP1","SNW1","TXNRD2","TXNRD3","TXNRD1","actr2","actr2b","C4A","CAPNS1",
				"CASP1","CDC42","CDK2","CHUK","COX7A3","CSF2RB","CSNK2A1","CYP11B2","CYP27C1","CYP2A6","CYP2A7",
				"CYP3A43","CYP4F8","DIO1","DIO2","DIO3","DNMT3a","DNMT3b","GBP1","GPX3","GPX2","GPX4","GPX1",
				//,"","","","","",
//				"","","","","","","","","","","","",
//				"","","","","","","","","","","","",
//				"","","","","","","","","","","","",
//				"","","","","","","","","","","","",
//				"","","","","","","","","","","","",
//				"","","","","","","","","","","","",
//				"","","","","","","","","","","","",
//				"","","","","","","","","","","","",
//				"","","","","","","","","","","","",
//				"","","","","","","","","","","","",
				};
		
		for(int i = 0; i<uppercase.length;i++){
		
			String inputID = uppercase[i];
			
			String outputID = inputID.substring(0,1)+inputID.substring(1,inputID.length()).toLowerCase();
			
			spf(inputID,outputID);
			
		}
		
	}


	
	private static void continueBuildmp(Vector<String> kv, String sid) {
		
		if(sid.contains("'")) return;
		
		for(int i = 0; i<kv.size();i++){
			
			if(kv.get(i).length()<3) return;
			
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
	
	private static void spmp(String spcf) throws IOException{
		
		BufferedReader br = new BufferedReader(new FileReader(spcf));
		
		String ln;
		
		String[] arr;
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split("=");
			
			spf(arr[0],arr[1]);
			
		}
		
		br.close();
		
		mapUpcase();

	}



	private static void spf(String iName, String oName) { //iName is input , oName is output standard name.
		
		Vector<String> v = new Vector<String>();
		
		v.add(iName);
		
		continueBuildmp(v,oName);
		
		
		
	}

	public HashMap<String, HashMap<String, Vector<String>>> groupmp(
			HashMap<String, HashMap<String, Vector<String>>> wprgro) {

		System.out.println("Start the normalization of Groups-Genes\n");
		
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
		
		System.out.println("Finished the normalization of Groups-Genes\n");
		
		return gro;
	}		
	
	
}


