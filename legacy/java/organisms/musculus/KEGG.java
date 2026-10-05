package musculus;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import tools.*;


public class KEGG{
    
    private static String ln,k,gk,gc,v,pathID,pathName;
    
    private static String[] arr;
    
	/* The HashMap is used to store the information of which genes are in certain pathways. With pathway name as 
	 * the key, gene names stored in a HashMap as a value.  
	 */
	public static HashMap<String,HashMap<String,String>> phge = new HashMap<String,HashMap<String,String>>();
	
	/* The HashMap is used to store the information of which gene pairs are in certain pathways. With pathway 
	 * name as the key, Value is a HashMap, gene pair names (A>B) as key stored , relationship as value in 
	 * this HashMap as a value.
	 */
	public static HashMap<String,HashMap<String,String>> phgp = new HashMap<String,HashMap<String,String>>();    
   
	
	/* The HashMap is used to store the information of which genes in groups are in certain pathways. With pathway name as the key, 
	 * GroupID as the Key of HashMap as the value of phgro, while gene names stored in a Vector as a value of Value HashMap.  
	 */
	public static HashMap<String,HashMap<String,Vector<String>>> phgro = new HashMap<String,HashMap<String,Vector<String>>>();
	
	/*
	 * This HashMap is used to map the Entrez ID into the Gene Name(HGNC Symbol and MGI Symbol)
	 */
	public static HashMap<String,Vector<String>> ecmp = new HashMap<String,Vector<String>>();
	
	public static void main(String args[]) throws IOException{
		
		extract("musculus");
		
	}
	
    public static void extract(String orgs)throws IOException{
        
    	String idkg  	= orgs+File.separator+"source"+File.separator+"KEGG"+File.separator;
        
    	String odkg  	= orgs+File.separator+"extraction"+File.separator+"KEGG"+File.separator;
        
    	//String kgmap 	= orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"uniprot-Mus";
    	
    	String kgmap 	= orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"MGIsymbo-Entreze";
        
    	String EtrMGI 	= orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"MGI_EntrezGene.rpt";
        //HashMap<String,Vector<String>> ecmp   = EntSybMap(kgmap);
    	String ncbi 	= orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"NCBImapping";    	    	
    	
    	String kompf 	= orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"mouseKO";
    	
    	EntrSybNCBI(ncbi);
    	
    	System.out.println("Size of mapping file after ncbi mapping file: "+ecmp.size());
    	
    	//test.MappingContainers(ecmp);
    	
//     	EntrMGISyb(EtrMGI);
//     	
//    	System.out.println("Size of mapping file after MGI mapping file: "+ecmp.size());
    	
    	EntSybMap(kgmap);    	
    	
    	System.out.println("Size of mapping file after biomart mapping file: "+ecmp.size());
    	
    	KoSybMap(kompf);
    	
    	System.out.println("Size of mapping file after ko mapping file: "+ecmp.size());
        
        iterExtract(idkg,odkg);        

    }


	
	private static void iterExtract(String idkg, String odkg) throws IOException {   	     	 
    	 
    	 File dir = new File(idkg);
    	 
         String[] filesNames = dir.list();
         
         for(int i = 0; i<filesNames.length;i++){
        	 
             String id = filesNames[i];
             //System.out.println(id.substring(0, id.length()-4));
             pathID= id.substring(0, id.length()-4);
             
             Vector<String[]> pairv = new Vector<String[]>();
             
             //HashMap<String,Vector<String>> ecmp = ECmap(kgmap);                    
             
             HashMap<String,Vector<String>> edgmap = EdgesMap(idkg+pathID+".xml", pairv,ecmp);   
             
             HashMap<String,Vector<String>> gpmap  = buildGroups(idkg+pathID+".xml"); 
             
             StoreOutputGenePairs(pairv, odkg +pathID+"genepair.txt",edgmap,gpmap,ecmp);
             
             //HashMap egmap = buildEdges(idkg+pathID+".xml", pairv,ecmp); 
             
             //writepair(pairv, odkg +pathID+"genepair.txt",egmap,gpmap,ecmp);
             
         }    	 
		
	}

	private static void KoSybMap(String kompf) throws IOException{
		// this method is used to include more information to the ecmp container, koid as key, GeneSymbol as value in Vector.
		
		String ln ;
		String[] arr;
		
		BufferedReader br = new BufferedReader(new FileReader(kompf));
		
		while((ln= br.readLine())!=null){
			
			arr = ln.split("\t");
			 
			if(arr.length<2||arr[0].length()<2||arr[1].length()<3||!arr[1].contains("K")){
				continue;
			}
			
			//System.out.println(ln);
			
			if(ecmp.containsKey(arr[0])){
				
				//System.out.println(arr[0]);
				
				Vector<String> sbv = ecmp.get(arr[0]);
				
				Vector<String> tv = new Vector<String>();
				
				//kv.add(arr[]);
				
				for(String sb:sbv){
					
					//System.out.println(sb);
					
					if(ecmp.containsKey(arr[1])){
						
						tv = ecmp.get(arr[1]);
						
						if(!tv.contains(sb)){
							
							tv.add(sb);
							
							ecmp.put(arr[1], tv);
							
						}
						
						
					}else{
						
						tv.add(sb);
						
						ecmp.put(arr[1], tv);
						
						//System.out.println(arr[1]+"\t"+sb);
						
					}
					
					
				}
				
				
				
			}
			
			
		}
		
		
		br.close();
		
	}

	

	private static void StoreOutputGenePairs(Vector<String[]> pairv,String dirOut, HashMap<String, Vector<String>> egmap,
			HashMap<String, Vector<String>> gpmap,HashMap<String, Vector<String>> ecmp) throws IOException{
		/*
		 *  Since the PPrel and ECrel have higher priority, therefore store them into the Pathway-GenePairs HashMap first.
		 *  Because StoreGenePairs() method implemented in a way that later input stuff will not not be stored if 
		 *  key already contained. 
		 */
		
        
        PrintWriter pw = new PrintWriter (new BufferedWriter(new FileWriter(dirOut)));
        
        Vector<String> va = new Vector<String>();
        
        Vector<String> vb = new Vector<String>();
        
        HashMap<String,Boolean> ck = new HashMap<String,Boolean>();
        
        String ita, itb;

        
        for(int z = 0;z<pairv.size();z++){
            
        	arr = pairv.get(z);
            
            if(egmap.containsKey(arr[0]) && egmap.containsKey(arr[1])){
                
            	ita = "";
            	itb = "";
            	                	
            	va = (Vector<String>) egmap.get(arr[0]);
                
                vb = (Vector<String>) egmap.get(arr[1]);    
                
                for(int r = 0; r < va.size(); r++){
                	
                	ita = va.get(r);
                	
                	for(int s = 0; s < vb.size();s++){                		                        
                        
                        itb = vb.get(s);                                            
                        
                        //WriteOutStoreInO(ita,itb,arr[2],pw,ecmp,ck);
                        
                        DirectWriteStore(ita,itb,arr[2],pathName,pw);
                		
                	}                	
                }                             
            }
        }
        
        //System.out.println("Debug: Size of gpmap"+ gpmap.size());
        
        
        if(gpmap.size()>0){
        	
        	HashMap<String,Vector<String>> vmp = new HashMap<String,Vector<String>>();
        	
            Set<String> gpset = gpmap.keySet();
            
            for(String gps : gpset){
            	                
            	ita = "";
            	
            	itb = "";
            	            	
            	if(gps.length()<1) continue;
            	
            	Vector<String> cv  = gpmap.get(gps);
            	
            	Vector<String> ncv = new Vector<String>();
            	
            	for(int i = 0; i< cv.size(); i++){
            		
            		va = (Vector<String>) egmap.get(cv.get(i));
            		
            		//System.out.println(pathName+"\t"+cv.get(i)+"\t"+gps);
            		//it shows at this step it is correct.
            		
            		for(String vas:va){
            			
            			//System.out.println(pathName+"\t"+vas+"\t"+gps);
            			//Also very precise here.
                    	
            			if(!ncv.contains(vas)){
                    		
                    		ncv.add(vas);
                    		
                    	}
            			
            		}
            		
            		for(int j = i+1;j<cv.size();j++){            		                    	
                        
                        vb = (Vector<String>) egmap.get(cv.get(j));                                                 
                        
                        //System.out.println(cv.get(i)+"--"+cv.get(j));
                        
                        for(int r = 0; r < va.size(); r++){
                        	
                        	ita = va.get(r);
                        	
//                        	if(!ncv.contains(ita)){
//                        		
//                        		ncv.add(ita);
//                        		
//                        	}
                        	
                        	for(int s = 0 ; s < vb.size();s++){                        	                                
                                
                                itb = vb.get(s);                                            
                                
                                //WriteOutStoreInO(ita,itb,"GPrel",pw,ecmp,ck);
                                
                                DirectWriteStore(ita,itb,"GPrel",pathName,pw);
                                
                                //System.out.println(ita+"\t"+itb+"\t"+"GPrel"+"\t"+pathName);
                                //Also correct here.
                        		
                        	}                	
                        }                         

            		}            		
            		
            	}
            	
            	if(!vmp.containsKey("Group"+gps)){
            		
            		vmp.put("Group"+gps,ncv);
            		
            	}else{
            		
            		System.out.println("Alert! It seems some group ids are the same!");
            		
            		ncv.addAll(vmp.get("Group"+gps));
            		
            		vmp.put("Group"+gps,ncv);
            		
            	}
            	
            }
        	
            phgro.put(pathName, vmp);
        }                
        pw.close();        
		
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
			
			if(ln.contains("Official Symbol:")){
				
				Pattern sbp = Pattern.compile("Official Symbol\\:(.*?)\\sand");
				
				Matcher sbm = sbp.matcher(ln);
				
				while(sbm.find()){
					
					symbo = sbm.group(1);
					
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
		
		/*
		System.out.println("The ncbi mapping file have just finish construction, size: "+ecmp.size()+
				"And how many lines of the mapping file: "+ indx);
		*/
	}
	

	private static void EntSybMap(String kgmap) throws IOException {
		// TODO Auto-generated method stub
		//HashMap<String, Vector<String>> esmp = new HashMap<String,Vector<String>>();								
		
		String[] arr;
		
		String ln;		
		
		BufferedReader br = new BufferedReader(new FileReader(kgmap));
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split("\t");
			
			if(arr.length<2)continue;
			
			if(arr[1].length()<2) continue;
			
			Vector<String> kv = new Vector<String>();
					
			if(arr[1].length()>1){
				
				kv.add(arr[1].trim());
			}
			
			
			if(arr[0].length()>1){
				
				if(ecmp.containsKey(arr[1])){
					/*
					vv = ecmp.get(kv.get(i));
					
					if(!vv.contains(sid)){
						
						vv.add(sid);
						
						ecmp.put(kv.get(i),vv);
					}
					*/
				}else{
					
					Vector<String> tvv = new Vector<String>();
					
					ecmp.put(arr[1],tvv);
					
				}
				
			}
						
		}
		
		
		br.close();		
		
		//return esmp;
	}
	
	
    private static void EntrMGISyb(String etrMGI) throws IOException{
		//This is to processing another entrez to MGI symbol mapping file, also store information into ecmp.
    	
		String[] arr;
		
		String ln;		
		
		BufferedReader br = new BufferedReader(new FileReader(etrMGI));
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split("\t");
			
			if(arr.length<9) continue;
			
			if(arr[1].length()<2) continue;
			
			if(arr[2].equalsIgnoreCase("W"))continue;
			
			Vector<String> kv = new Vector<String>();
					
			if(arr[1].length()>1){
				
				kv.add(arr[1].trim());
				
			}
			
			
			if(arr[0].length()>1){
				
				continueBuildmp(kv, arr[0].trim());
				
			}
						
		}
		
		
		br.close();	
		
	}
    
	
	
	
	
	private static void continueBuildmp(Vector<String> kv, String sid) {
		
		//if(sid.contains("'")) return;
		
		for(int i = 0; i<kv.size();i++){
			
			if(kv.get(i).length()<1) return;
			
			Vector<String> vv = new Vector<String>();
			
			if(ecmp.containsKey(kv.get(i))){
				/*
				vv = ecmp.get(kv.get(i));
				
				if(!vv.contains(sid)){
					
					vv.add(sid);
					
					ecmp.put(kv.get(i),vv);
				}
				*/
			}else{
				
				vv.add(sid);
				
				ecmp.put(kv.get(i),vv);
				
			}
			
		}		
		
	}

		
	//This is the newly implemented function that based only on the KEGG entry to get genes rather than KEGG IDs.
	
	private static HashMap<String,Vector<String>> EdgesMap(String dir, Vector<String[]> pairv, HashMap<String,Vector<String>> ecmp) throws IOException {
        
		HashMap<String, Vector<String>> edgmap = new HashMap<String, Vector<String>>(); 
        		
		BufferedReader br = new BufferedReader(new FileReader(dir));
        
		String key;
        
		k = "";        
        
        while((ln = br.readLine())!=null){
            
            if(ln.contains("title=")){
                
            	Pattern pn = Pattern.compile("title=\"(.*)\"");
                
            	Matcher mn = pn.matcher(ln);
                
            	while(mn.find()){
                
            		pathName = mn.group(1);
            		
                }
            }
            
            Vector<String> val = new Vector<String>();
            
            if(ln.contains(" <entry id=")&&(ln.contains("type=\"gene\""))){
        
            	Pattern pe = Pattern.compile("<entry id=\"(\\d*)\" name=");
                
            	Matcher m = pe.matcher(ln);
                
            	while(m.find()){
            		
                    k = m.group(1);
                    
                }
            	
            	Vector<String> nv = new Vector<String>();
            	
            	Pattern pn = Pattern.compile("mmu:(\\d*)");
                
            	Matcher mn = pn.matcher(ln);
            	
                //int i = 1;
                
            	while(mn.find()){
            		
                    String entrz = mn.group(1).trim();
                    
                    //i++;  
            		//String geneName;      
            		
            		String pathway = pathName;
            		
            		Vector<String> tmkv = new Vector<String>();
                	
                	if(ecmp.containsKey(entrz)){     
                                   
                		tmkv = ecmp.get(entrz);
                		
                        for(String tmk : tmkv){
                        	
                            if(!nv.contains(tmk)){
                            	
                            	nv.add(tmk);
                            	                                   	
                            }
                            
                            StoreGenes(tmk,pathway);
                        	
                        }
                		
                	
                	} 
                	/*
                	else{
                		
                		//If don't want too many Entrez ID can't be mapped just block the following part. 
                		System.out.println("There are some Entrez ID haven't been mapped: "+entrz);                		                		
                		
                		StoreGenes(entrz,pathway);
                		 
                		 
                		if(!nv.contains(entrz)){
                            	
                            	nv.add(entrz);
                            	                                   	
                        }
                		
                		
                	} 
                    */                                                                                                                     
                }          
            	
            	Vector<String> tvv = new Vector<String>();
            	
            	//if(k.length()<1) continue;
            	
            	//System.out.println("Debug: the number of corresponding Entrez ID to a KEGG entry: "+nv.size());
            	
            	if(edgmap.containsKey(k)){
            		
            		tvv = edgmap.get(k);            		
            		
            		for(String ns: nv){
            			
            			if(!tvv.contains(ns)){
            				
            				tvv.add(ns);
            				
            			}
            		}
            		
            		edgmap.put(k, tvv);
            		
            		
            	}else{
            		
            		edgmap.put(k, nv);
            		
            	}
            	

                
                k  = ""; 
                
                nv = new Vector<String>();
                
            }
            
            
            if(ln.contains(" <entry id=")&&(ln.contains("type=\"ortholog\""))){
                
            	Pattern pe = Pattern.compile("<entry id=\"(\\d*)\" name=");
                
            	Matcher m = pe.matcher(ln);
                
            	while(m.find()){
            		
                    k = m.group(1);
                    
                }
            	
            	Vector<String> nv = new Vector<String>();
            	
            	Pattern pn = Pattern.compile("ko:(K\\d*)");
                
            	Matcher mn = pn.matcher(ln);
            	
                //int i = 1;
                
            	while(mn.find()){
            		
                    String entrz = mn.group(1).trim();
                    
                    //System.out.println(entrz);
                    
                    //i++;  
            		//String geneName;      
            		
            		String pathway = pathName;
            		
            		Vector<String> tmkv = new Vector<String>();
                	
                	if(ecmp.containsKey(entrz)){     
                                   
                		//System.out.println(entrz);
                		
                		tmkv = ecmp.get(entrz);
                		
                        for(String tmk : tmkv){
                        	
                            if(!nv.contains(tmk)){
                            	
                            	//System.out.println(tmk);
                            	
                            	nv.add(tmk);
                            	                                   	
                            }
                            
                            StoreGenes(tmk,pathway);
                        	
                        }
                		
                	
                	}                                                                                                                      
                }          
            	
            	Vector<String> tvv = new Vector<String>();
            	
            	//if(k.length()<1) continue;
            	
            	//System.out.println("Debug: the number of corresponding Entrez ID to a KEGG entry: "+nv.size());
            	
            	if(edgmap.containsKey(k)){
            		
            		tvv = edgmap.get(k);            		
            		
            		for(String ns: nv){
            			
            			if(!tvv.contains(ns)){
            				
            				tvv.add(ns);
            				
            			}
            		}
            		
            		edgmap.put(k, tvv);
            		
            		
            	}else{
            		
            		edgmap.put(k, nv);
            		
            	}
            	

                
                k  = ""; 
                
                nv = new Vector<String>();
                
            }
   
            
            String[] tmp = new String[3];            
            
            if(ln.contains("<relation entry1=")){
                //System.out.println("debuging1");
            
            	Pattern pp = Pattern.compile("<relation entry1=\"(\\d*)\" entry2=\"(\\d*)\" type=\"(\\w*)\">");
                
            	Matcher mp = pp.matcher(ln);
                
            	while(mp.find()){
            		
                    tmp[0] = mp.group(1);

                    tmp[1] = mp.group(2);
                    
                    tmp[2] = mp.group(3);
                    
                    pairv.add(tmp);
                }
                
            }
        }
        
        br.close();
        
        return edgmap;
    }


    
	private static HashMap<String, Vector<String>> buildGroups(String dir) throws IOException{
		// This method is specifically extract the group information.
		
		HashMap<String,Vector<String>> gpmp = new HashMap<String,Vector<String>>();
		
		Vector<String> vv = new Vector<String>();
		
		BufferedReader br = new BufferedReader(new FileReader(dir));
        
		String key;
        k = "";
        gk = "";
        
        while((ln = br.readLine())!=null){         
            
            
            if(ln.contains(" <entry id=")&&(ln.contains("type=\"group\""))){
                
            	Pattern pe = Pattern.compile("<entry id=\"(\\d*)\"");
                
            	Matcher mgp = pe.matcher(ln);
                
            	while(mgp.find()){
            		
                    gk = mgp.group(1);
                    
                    //System.out.println("Debug of missing group key"+gk);
                    
                }                                
            }
            
            if(ln.contains("<component id=")){
                
            	Pattern ce = Pattern.compile("<component id=\"(\\d*)\"");
                
            	Matcher mce = ce.matcher(ln);
                
            	while(mce.find()){
            		
                    gc = mce.group(1);                                        
                    
                    if(!vv.contains(gc) && gc!=""){
          	
                    	Vector<String> tv = new Vector<String>();
                    	
                    	//System.out.println("Debug of missing group component"+gc);
                    	
                    	if(gk.length()>0 && gk!=""){
                    		
                        	if(gpmp.containsKey(gk)){
                        		
                        		tv = gpmp.get(gk);
                        		
                        		if(!tv.contains(gc)){
                        			
                        			tv.add(gc);
                        			
                        			gpmp.put(gk, tv);
                        			
                        		}
                        		
                        	}else{
                        		
                        		tv.add(gc);
                        		
                        		gpmp.put(gk, tv);
                        		
                        	}	
                    		
                    		
                    	}
                    	
                    	tv = new Vector<String>();
                    }                                         
                    
                    
                    gc = "";
                    
                    vv = new Vector<String>();
                }                                
            }
            

            
            if(ln.contains("</entry>")){             	                            	                
                
                k  = "";  
                
                gk = "";
                
                gc = "";
                
                vv = new Vector<String>();
                
            }
           
        }
        
        br.close();
        
		return gpmp;
	}

	private static HashMap buildegmap(String dir, Vector<String[]> pairv, HashMap<String,Vector<String>> ecmp) throws IOException {
        
		HashMap<String, Vector> egmap = new HashMap<String, Vector>(); 
        		
		BufferedReader br = new BufferedReader(new FileReader(dir));
        
		String key;
        
		k = "";
        
        //gk = "";
        
        while((ln = br.readLine())!=null){
            
            if(ln.contains("title=")){
                
            	Pattern pn = Pattern.compile("title=\"(.*)\"");
                
            	Matcher mn = pn.matcher(ln);
                
            	while(mn.find()){
                
            		pathName = mn.group(1);
            		
                }
            }
            
            Vector<String> val = new Vector<String>();
            
            if(ln.contains(" <entry id=")&&(ln.contains("type=\"gene\""))){
        
            	Pattern pe = Pattern.compile("<entry id=\"(\\d*)\" name=");
                
            	Matcher m = pe.matcher(ln);
                
            	while(m.find()){
            		
                    k = m.group(1);
                    
                }
                //System.out.println(k);                 
            }
            
            if(ln.contains("<graphics name=")){
                
            	Pattern pg = Pattern.compile("<graphics name=\"(.*)\" fgcolor=");
                
            	Matcher m = pg.matcher(ln);
                
            	while(m.find()){
                
            		v = m.group(1);                    
                    //System.out.println(v);
                }                
            }
            
            if(ln.contains("</entry>")){                
                arr = v.split(", ");
                /*
                for(int i = 0;i<arr.length;i++){
                    if(v.contains("...") && i == arr.length-1){
                        arr[i] = arr[i].substring(0, arr[i].length()-3);
                        
                    }
                    val.add(arr[i]);
                }        
                 * This is a bug caused by the miss understanding of KGML
                 */
                if(arr[0].contains("...")){
                    arr[0] = arr[0].substring(0, arr[0].length()-3);   
                }
             
                val.add(arr[0]);                   
                
                //System.out.println(val.size());
                if(k!="" && val.size()>0 && !egmap.containsKey(k)){
                    
                	egmap.put(k, val);
                    
                    val = new Vector<String>();
                    
                }
                
                if(k!=""){
                	
            		String geneName;      
            		
            		String pathway = pathName;
                	
                	if(ecmp.containsKey(arr[0])){                	
                		
                		Vector<String> genamV =ecmp.get(arr[0]);
                		
                		for (Object gno : genamV){
                			
                			geneName = gno.toString();
                			
                			StoreGenes(geneName,pathway);
                			
                		}                		                	
                	
                	} else{
                		
                		geneName = arr[0];
                		
                		StoreGenes(geneName,pathway);
                	}                	               	
                	
                }
                
                k  = "";  
                //gk = "";
                
            }
   
            
            String[] tmp = new String[3];            
            
            if(ln.contains("<relation entry1=")){
                //System.out.println("debuging1");
            
            	Pattern pp = Pattern.compile("<relation entry1=\"(\\d*)\" entry2=\"(\\d*)\" type=\"(\\w*)\">");
                
            	Matcher mp = pp.matcher(ln);
                
            	while(mp.find()){
            		
                    //System.out.println("debuging2");
                    tmp[0] = mp.group(1);
                    //System.out.println(tmp[0]);
                    
                    tmp[1] = mp.group(2);
                    
                    tmp[2] = mp.group(3);
                    
                    pairv.add(tmp);
                }
                
            }
        }
        
        br.close();
        
        return egmap;
    }

	private static HashMap buildEdges(String dir, Vector<String[]> pairv, HashMap<String,Vector<String>> ecmp) throws IOException {
        
		HashMap<String, String> egmap = new HashMap<String, String>(); 
        		
		BufferedReader br = new BufferedReader(new FileReader(dir));
        
		String key;
        k = "";
        
        while((ln = br.readLine())!=null){
            
            if(ln.contains("title=")){
                
            	Pattern pn = Pattern.compile("title=\"(.*)\"");
                
            	Matcher mn = pn.matcher(ln);
                
            	while(mn.find()){
                
            		pathName = mn.group(1);
            		
                }
            }
            
            //Vector<String> val = new Vector<String>();
            
            if(ln.contains(" <entry id=")&&(ln.contains("type=\"gene\""))){
        
            	Pattern pe = Pattern.compile("<entry id=\"(\\d*)\" name=");
                
            	Matcher m = pe.matcher(ln);
                
            	while(m.find()){
            		
                    k = m.group(1);
                    
                }
                //System.out.println(k);                 
            }
            
            if(ln.contains("<graphics name=")){
                
            	Pattern pg = Pattern.compile("<graphics name=\"(.*)\" fgcolor=");
                
            	Matcher m = pg.matcher(ln);
                
            	while(m.find()){
                
            		v = m.group(1);                    
                    //System.out.println(v);
                }                
            }
            
            if(ln.contains("</entry>")){                
                arr = v.split(", ");
                /*
                for(int i = 0;i<arr.length;i++){
                    if(v.contains("...")&& i == arr.length-1){
                        arr[i] = arr[i].substring(0, arr[i].length()-3);
                        
                    }
                    val.add(arr[i]);
                }        
                 * This is a bug caused by the miss understanding of KGML
                 */
                if(arr[0].contains("...")){
                    arr[0] = arr[0].substring(0, arr[0].length()-3);   
                }
             
                egmap.put(k, arr[0]);
                
                if(k!=""){
                	
            		String geneName;      
            		
            		String pathway = pathName;
                	
                	if(ecmp.containsKey(arr[0])){                	
                		
                		Vector<String> genamV =ecmp.get(arr[0]);
                		
                		for (Object gno : genamV){
                			
                			geneName = gno.toString();
                			
                			StoreGenes(geneName,pathway);
                			
                		}                		                	
                	
                	} else{
                		
                		geneName = arr[0];
                		
                		StoreGenes(geneName,pathway);
                	}                	               	
                	
                }
                
                k="";                                                                        
                
            }
   
            
            String[] tmp = new String[3];            
            
            if(ln.contains("<relation entry1=")){
                //System.out.println("debuging1");
            
            	Pattern pp = Pattern.compile("<relation entry1=\"(\\d*)\" entry2=\"(\\d*)\" type=\"(\\w*)\">");
                
            	Matcher mp = pp.matcher(ln);
                
            	while(mp.find()){
            		
                    //System.out.println("debuging2");
                    tmp[0] = mp.group(1);
                    //System.out.println(tmp[0]);
                    
                    tmp[1] = mp.group(2);
                    
                    tmp[2] = mp.group(3);
                    
                    pairv.add(tmp);
                }
                
            }
        }
        
        br.close();
        
        return egmap;
    }
	
    private static void StoreGenes(String geneName,String pathway) {
        
    	HashMap<String, String> ges = new HashMap<String,String>();  
        
        //if (arr[0].length()<1) continue;
        
        if (phge.containsKey(pathway)){
        	
        	ges = phge.get(pathway);
            
        	if (!ges.containsKey(geneName)){
            	
        		ges.put(geneName,"KEGG");
        		
            	phge.put(pathway, ges);        		
            	
            }            	        	
        	
        } else {
        	
        	ges.put(geneName, "KEGG");
        	
        	phge.put(pathway, ges);
        	
        }  
		
	}

	private static void writepair(Vector<String[]> pairv, String dirOut, HashMap egmap,HashMap<String,Vector<String>> gpmap,HashMap<String,Vector<String>> ecmp) throws IOException{
            
            PrintWriter pw = new PrintWriter (new BufferedWriter(new FileWriter(dirOut)));
            
            Vector<String> va = new Vector<String>();
            
            Vector<String> vb = new Vector<String>();
            
            HashMap<String,Boolean> ck = new HashMap<String,Boolean>();
            
            String ita, itb;

            
            for(int z = 0;z<pairv.size();z++){
                
            	arr = pairv.get(z);
                
                if(egmap.containsKey(arr[0])&&egmap.containsKey(arr[1])){
                    
                	ita = "";
                	itb = "";
                	                	
                	va = (Vector<String>) egmap.get(arr[0]);
                    
                    vb = (Vector<String>) egmap.get(arr[1]);                    
                    
                    ita = va.get(0);
                    itb = vb.get(0);
                                        
                    /*
                    ita = (String) egmap.get(arr[0]);
                    
                    itb = (String) egmap.get(arr[1]);
                    */
                    
                    //WriteOutStoreInO(ita,itb,arr[2],pw,ecmp,ck);
                    
                    DirectWriteStore(ita,itb,arr[2],pathName,pw);
                    
                }
            }
            
            //System.out.println("Debug: Size of gpmap"+ gpmap.size());
            
            if(gpmap.size()>0){
            	
                Set<String> gpset = gpmap.keySet();
                
                for(String gps : gpset){
                	                
                	ita = "";
                	
                	itb = "";
                	
                	
                	if(gps.length()<1) continue;
                	
                	Vector<String> cv = gpmap.get(gps);
                	
                	for(int i = 0; i< cv.size(); i++){
                		
                		for(int j = i+1;j<cv.size();j++){
                		
                        	va = (Vector<String>) egmap.get(cv.get(i));
                            
                            vb = (Vector<String>) egmap.get(cv.get(j));                    
                            
                            //System.out.println(cv.get(i)+"--"+cv.get(j));
                            
                            ita = va.get(0);
                            
                            itb = vb.get(0);                			                			                		
                			
                			//WriteOutStoreInO(ita,itb,"GPrel",pw,ecmp,ck);
                			
                			DirectWriteStore(ita,itb,"GPrel",pathName,pw);
                		}
                		
                		
                	}
                	
                }
            	
            }
            

            
            pw.close();        
    }

    private static void DirectWriteStore(String ita, String itb, String rel,
			String pwy, PrintWriter pw) throws IOException{
		// TODO Auto-generated method stub
    	
    	StoreGenePairs(ita, itb, rel,pwy);
    	
    	pw.println(ita+"\t"+itb+"\t"+rel+"\t"+pwy);
		
	}

	private static void WriteOutStoreInO(String ita, String itb, String rel,PrintWriter pw, HashMap<String, Vector<String>> ecmp, HashMap<String, Boolean> ck) {
		// TODO Auto-generated method stub
    	
	  if(!ecmp.containsKey(ita)&&!ecmp.containsKey(itb)){
	      
	  if(!ck.containsKey(ita+"\t"+itb+"\t"+arr[2]+pathName)&&!ck.containsKey(itb+"\t"+ita+"\t"+arr[2]+pathName)){
	          pw.println(ita+"\t"+itb+"\t"+arr[2]+"\t"+pathName+"\t"+pathID);
	          ck.put(ita+"\t"+itb+"\t"+arr[2], Boolean.TRUE);                                                        
	      }
	  	StoreGenePairs(ita, itb, rel,pathName);
	      
	  }else if(ecmp.containsKey(ita)&&!ecmp.containsKey(itb)){
	      
	  	Vector<String> av = ecmp.get(ita);
	      
	  	for(int a = 0;a<av.size();a++){
	          String itag = av.get(a);
	          if(!ck.containsKey(itag+"\t"+itb+"\t"+rel+pathName)&&!ck.containsKey(itb+"\t"+itag+"\t"+rel+pathName)){
	              pw.println(itag+"\t"+itb+"\t"+rel+"\t"+pathName+"\t"+pathID);
	              ck.put(itag+"\t"+itb+"\t"+rel+pathName, Boolean.TRUE);                                                               
	          }
	          StoreGenePairs(itag, itb, rel,pathName);
	      }
	  	
	  }else if(ecmp.containsKey(itb)&&!ecmp.containsKey(ita)){
	      Vector<String> bv = ecmp.get(itb);
	      for(int b = 0;b<bv.size();b++){
	          String itbg = bv.get(b);
	          if(!ck.containsKey(itbg+"\t"+ita+"\t"+rel+pathName)&&!ck.containsKey(ita+"\t"+itbg+"\t"+rel+pathName)){
	              pw.println(itbg+"\t"+itb+"\t"+rel+"\t"+pathName+"\t"+pathID);
	              ck.put(itbg+"\t"+itb+"\t"+rel+pathName, Boolean.TRUE);                                                                
	          }
	          StoreGenePairs(itbg, itb, rel,pathName);
	      }
	  }else{
	      Vector<String> av = ecmp.get(ita);
	      Vector<String> bv = ecmp.get(itb);
	      for(int a = 0;a<av.size();a++){
	          for(int b = 0;b<bv.size();b++){
	              String itag = av.get(a);
	              String itbg = bv.get(b);
	              if(!ck.containsKey(itbg+"\t"+itag+"\t"+rel+pathName)&&!ck.containsKey(itag+"\t"+itbg+"\t"+rel+pathName)){
	                  pw.println(itag+"\t"+itbg+"\t"+rel+"\t"+pathName+"\t"+pathID);
	                  ck.put(itag+"\t"+itbg+"\t"+rel+pathName, Boolean.TRUE);                                                                        
	              }
	              StoreGenePairs(itag, itbg, rel,pathName);
	          }
	      }
	
	  }
		
	}
    
	private static void StoreGenePairs(String geneA, String geneB, String rel,String pathway) {
		
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

	private static HashMap<String, Vector<String>> ECmap(String idirmp) throws IOException {
        
		String ln;
        
        String[] arr;

        HashMap<String,Vector<String>> ecmp = new HashMap<String,Vector<String>>();
        
        BufferedReader br = new BufferedReader(new FileReader(idirmp));
        
        while((ln = br.readLine())!=null){

            if(ln.contains("EC")){
                
            	arr = ln.split("\t");
                
                String[] ge = arr[4].split(" ");

                Pattern ep = Pattern.compile("EC (\\d{1,4}\\.\\d{0,4}\\.\\d{0,4}\\.[\\d|-]{0,4})");
                
                Matcher em = ep.matcher(ln);
                
                while(em.find()){
                    //System.out.println(em.group(1));
                    
                	Vector<String> vec = new Vector<String>();
                    
                    if(!ecmp.containsKey(em.group(1))){
                        
                    	vec.add(ge[0]);
                        
                        ecmp.put(em.group(1),vec);
                        
                        ecmp.put(em.group(),vec);
                    
                    }else if(ecmp.containsKey(em.group(1))){
                        
                    	vec= ecmp.get(em.group(1));
                        
                        if(!vec.contains(ge[0])){
                        
                            vec.add(ge[0]);
                            
                            ecmp.put(em.group(1),vec);
                            
                            ecmp.put(em.group(),vec);
                        	
                        }
                    }
                }
            }
        }
        //System.out.println("The size of EC-Gene Name mapping: "+ecmp.size());
        br.close();
        return ecmp;
    }
}
