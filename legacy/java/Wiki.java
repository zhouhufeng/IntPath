/* Title: WikiPathway parser
 * Function: Java Program to Extrace Information from GPML 
 * 
 * Author: Hufeng Zhou
 * Time: July 3rd 2021
 * Version: v2
 */

package sapiens;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class Wiki{

    private static String ln,k,v,pathID,pathName;
    
    private static String[] arr;	
    
    public static HashMap<String,String> wikimp = new HashMap<String,String>();
    		
	/* The HashMap is used to store the information of which genes are in certain pathways. With pathway name as 
	 * the key, gene names stored in a HashMap as a value.  
	 */
	public static HashMap<String,HashMap<String,String>> phge = new HashMap<String,HashMap<String,String>>(10000);
	
	
	/* The HashMap is used to store the information of which gene pairs are in certain pathways. With pathway 
	 * name as the key, Value is a HashMap, gene pair names (A>B) as key stored , relationship as value in 
	 * this HashMap as a value.
	 */
	public static HashMap<String,HashMap<String,String>> phgp = new HashMap<String,HashMap<String,String>>(10000);    
	
	
	/* The HashMap is used to store the information of which genes in groups are in certain pathways. With pathway name as the key, 
	 * GroupID as the Key of HashMap as the value of phgro, while gene names stored in a Vector as a value of Value HashMap.  
	 */
	public static HashMap<String,HashMap<String,Vector<String>>> phgro = new HashMap<String,HashMap<String,Vector<String>>>();

	
	
	public static void extract(String orgs)throws IOException{
		
    	String idwp  = orgs+File.separator+"source"+File.separator+"WikiPathways"+File.separator;
    	
        String odwp  = orgs+File.separator+"extraction"+File.separator+"WikiPathways"+File.separator;
        
//        String idwp  = orgs+File.separator+"source"+File.separator+"Oldwiki"+File.separator;
//    	
//        String odwp  = orgs+File.separator+"extraction"+File.separator+"Oldwiki"+File.separator;
		
        String wpmap = orgs+File.separator+"source"+File.separator+"mapping"+File.separator+"uniprot-Mus";
        
        iterExtract(idwp,odwp,wpmap);
		
	}


	private static void iterExtract(String idwp, String odwp, String wpmap) throws IOException{		
		
        File folder = new File(idwp);
        
        String[] fileNames = folder.list();
        //System.out.println(fileNames.length);
        
        for(int i = 0; i<fileNames.length; i++){
            
            String id = fileNames[i];
            
            //System.out.println(id.substring(0, id.length()-5));
            pathID = id.substring(0, id.length()-5);
            
            ///*
            //This Vector contains String[] which stores the element pairs of the graph ID of a line entry.
            Vector<ArrayList<String>> linepairv = new Vector<ArrayList<String>>();
            
            //Construct the HashMap nodemp, using the graph ID as Key, the gene name and entreze ID put in a vector as value.
            HashMap<String,ArrayList> gragenmp = new HashMap<String,ArrayList>();
            
            //use the group ID as key while the Vector contains other gene name and Entrze ID as the value
            HashMap<String,Vector> groupmp = new HashMap<String,Vector>();
            
            //System.out.println(pathName);
            
            pathName = "";
            
            buildegmap(idwp+pathID+".gpml", linepairv, gragenmp, groupmp);
            //System.out.println(gragenmp.size());
            //System.out.println(linepairv.size());                       
            
            writepair(linepairv, gragenmp, groupmp, odwp+pathID+"genepair.txt");
            //*/
            
            /*
            pathName = "";
            
            Vector<String[]> lnpv = new Vector<String[]>();
            
            HashMap<String,HashMap<String,String>> gropmp = new HashMap<String,HashMap<String,String>>();
            
            HashMap<String,String[]> gaphmp = new HashMap<String,String[]>();
            
            ConstructEdges(idwp+pathID+".gpml",lnpv,gaphmp,gropmp);
            
            OutPutEdges(lnpv,gaphmp,gropmp,odwp+pathID+"genepair.txt");
            */
            
        }
	}


    private static void OutPutEdges(Vector<String[]> lnpv, HashMap<String, String[]> gaphmp,
			HashMap<String, HashMap<String, String>> groupmp, String odir) throws IOException{

        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(odir)));
        
        HashMap<String,Boolean> ckamp = new HashMap<String,Boolean>();
        
        HashMap<String,Boolean> ckbmp = new HashMap<String,Boolean>();

        for(int i = 0; i<lnpv.size();i++){
            
        	String[] gaid = lnpv.get(i);
            
        	String na = gaid[0];
            
        	String nb = gaid[1];
            
        	if(gaphmp.containsKey(na) && gaphmp.containsKey(nb)){
                
        		String[] genA = gaphmp.get(na);
                
        		String[] genB = gaphmp.get(nb);
        		
//        		if(genA[1]!=null && genA[1].length()>1){
//        			
//        			wikimp.put(genA[0],genA[1]);
//        		}
//        		
//        		if(genB[1]!=null && genB[1].length()>1){
//        			
//        			wikimp.put(genB[0],genB[1]);
//        		}
        		
                //System.out.println(genA.get(0)+"\t"+"|"+"\t"+genB.get(0) +"\t"+"|"+"\t"+"graphRel"+"\t"+"|"+"\t"+pathName+"\t"+"|"+"\t"+genA.get(1)+"\t"+"|"+"\t"+genB.get(1));
                
        		if(genA[0].length()>0 && genB[0].length()>0){//genA.get(0).length()<10 && genB.get(0).length()<10 && 
                    
        			if(!ckamp.containsKey(genA[0]+"\t"+genB[0])&&!ckamp.containsKey(genB[0]+"\t"+genA[0])){
                                				                       
        				ckamp.put(genA[0]+"\t"+genB[0], Boolean.TRUE);
        				
        				String gepreA = genA[0];
        				
        				String geA = gepreA.replaceAll("\\?", "").replaceAll("`","");//.replaceAll("@", "")
        				
        				String gepreB = genB[0];
        				
        				String geB = gepreB.replaceAll("\\?", "").replaceAll("`","");//.replaceAll("@", "")        				
        				        				
        				pw.println(geA+"\t"+geB +"\t"+"graphRel"+"\t"+pathName+"\t");
                        
        				StoreGenePairs(geA, geB, "graphRel", pathName);
        				
                        //System.out.println("hello>"+genA.get(0)+"\t"+genB.get(0));
                    }
                }
            }
        }
        
        
        Set<String> groupsets = groupmp.keySet();  
        
        System.out.println("The size of groups in each pathway: "+groupmp.size());
        
        for(String id:groupsets){            
        	
        	HashMap<String,String> egp = groupmp.get(id);
        	
        	Set<String> egpks = egp.keySet();
        	
        	Vector<String> env = new Vector<String>();
        	
        	for(String pk: egpks){
        		
        		if(pk.length()>1)
        			
        		env.add(pk);
        		
        	}
        	
        	System.out.println("The num of genes in each group: "+env.size());
        	for(int i = 0; i<env.size();i++){
                
        		for(int j = i+1; j<env.size();j++){ 
                    
        			String genea = (String) env.get(i);
                    String gena = genea.replaceAll("\\?", "").replaceAll("`","");//replaceAll(" ", "").toUpperCase()
        			
                    String geneb = (String) env.get(j); 
        			String genb = geneb.replaceAll("\\?", "").replaceAll("`","");//replaceAll(" ", "").toUpperCase()
                    
        			//System.out.println(gena+"\t"+"|"+"\t"+genb +"\t"+"|"+"\t"+"groupRel"+"\t"+"|"+"\t"+pathName);
                    
        			if(gena.length()<10&&genb.length()<10&&!gena.equalsIgnoreCase(genb)){
                        
        				if(!ckbmp.containsKey(gena+"\t"+genb)&&!ckbmp.containsKey(genb+"\t"+gena)){
                            
        					ckbmp.put(gena+"\t"+genb, Boolean.TRUE);
                            
        					pw.println(gena+"\t"+genb +"\t"+"groupRel"+"\t"+pathName);
                            
        					StoreGenePairs(gena, genb, "groupRel", pathName);
        					
                        }                       
                    }
                }
            }
            
        }
        
        pw.close();
		
	}


	private static void ConstructEdges(String dir, Vector<String[]> lnpv,HashMap<String, String[]> gaphmp,
			HashMap<String, HashMap<String, String>> groupmp) throws IOException {
		// TODO Auto-generated method stub
        BufferedReader br = new BufferedReader(new FileReader(dir));

        k = null;

        String gpk = "";
        
        int lidx = 0;
        
        String[] ge = new String[2];
        
        String[] lids = new String[2];               
        
        String namefg = "";

        while((ln = br.readLine())!=null){            
            //String[] gene = new String[2];

            if(ln.contains("<Pathway")){
            	
                Pattern pname = Pattern.compile("Name=\"(.*?)\"");
                
                Matcher mname = pname.matcher(ln);
                
                if(mname.find()){
                
                	pathName = mname.group(1);
                	
                }
            }
            

            
            if(ln.contains(" <DataNode ")){         
            	
            	HashMap<String,String> gpmpe = new HashMap<String,String>();
            	
            	//Vector<String> gphv = new Vector<String>();
            	
                Pattern genename = Pattern.compile("TextLabel=\"(.*?)\".*Type=\"GeneProduct\"");//

                Matcher mg = genename.matcher(ln);
                
                String genam = null;
                                
                
                if(mg.find()){
                	
                	String genamP = mg.group(1);
                	
                	genam = genamP.replaceAll("\\?", "").replaceAll("`","").replaceAll("@", "");//.toUpperCase()                	
                	
                	//gphv.add(genam);

                    //gpv.add(genam);
                    
                    ge[0]=genam;
                    
                    namefg = genam;
                    
                    gpmpe.put(genam, null);

                    StoreGenes(genam, pathName);

                }

                Pattern graphID = Pattern.compile("GraphId=\"(.*?)\"");
                
                Matcher mgra = graphID.matcher(ln);
                
                while(mgra.find()){
                	
                    k = mgra.group(1);

                }
                
                

                Pattern groupID = Pattern.compile("GroupRef=\"(.*?)\"");
                
                Matcher mgp = groupID.matcher(ln);

                if(mgp.find()){
                	
                    gpk = mgp.group(1);

                    //System.out.println(gpk);
                    if(!groupmp.containsKey(gpk)){
                    	
                    	//tmp.put(namefg, "null");
                    	
                        groupmp.put(gpk,gpmpe);
                        
                        //gpv = new Vector<String>();
                        
                    }else{
                    	HashMap<String,String> tmp = groupmp.get(gpk);
                        //if(mg.find()){
                    	
                    	if(genam!=null){

                            //System.out.println(mg.group(1));                            
                            
                            if(!tmp.containsKey(genam)){
                            	
                            	tmp.put(genam,null);
                            	
                            	genam = null;
                            	
                            	groupmp.put(gpk,tmp);
                            	
                            }
                                                                                    
                            //gpmpe = new HashMap<String,String>();

                    	}
                    }
                    
                    gpmpe = new HashMap<String,String>();
                    
                }
            }

            if(ln.contains("<Xref Database=")){
                
            	Pattern geneid = Pattern.compile("<Xref Database=\"(.*)\" ID=\"(\\w*)\" />");
                
                Matcher mgid = geneid.matcher(ln);
                
                while(mgid.find()){
                    //gv.set(1,mgid.group(1)+"\t"+mgid.group(2));
                	
                	ge[1] = mgid.group(2);
                	
                }                
                
            }

            if(ln.contains("</DataNode>")){
                
            	if(k!=null && (ge[0]!=null) && !gaphmp.containsKey(k)){
                    
            		gaphmp.put(k, ge);

                }                
            	
            	ge = new String[2];
            	
            }

            if(ln.contains("<Point")){
                //Pattern nodp = Pattern.compile("GraphRef=\"(\\w*)\"");
                Pattern nodp = Pattern.compile("GraphRef=\"(.*?)\"");
                
                Matcher nodm = nodp.matcher(ln);                                
                
                while(nodm.find()){
                	
                	lids[lidx] = nodm.group(1);

                	lidx++;                	
                	
                }
            }

            if(ln.contains("</Line>")){
            	
                if(lidx==2){
                	
                	lnpv.add(lids);
                	
                }
                
                lidx = 0;
                
                lids = new String[2];
                
            }

        }
        br.close();
    	
    	
	}


	private static void buildegmap(String dir, Vector<ArrayList<String>> linepairv, HashMap<String, ArrayList> gragenmp,
            HashMap<String, Vector> groupmp) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(dir));

        k = null;

        ArrayList<String> gv = new ArrayList<String>();
        gv.add(null); gv.add(null);

        ArrayList<String> nodv = new ArrayList<String>();

        String gpk = "";

        while((ln = br.readLine())!=null){            
            //String[] gene = new String[2];

            if(ln.contains("<Pathway")){
            	
                Pattern pname = Pattern.compile("Name=\"(.*?)\"");
                
                Matcher mname = pname.matcher(ln);
                
                if(mname.find()){
                
                	pathName = mname.group(1);
                    //System.out.println(pathName);
                	
                }
                //System.out.println(ln);
            }

            if(ln.contains(" <DataNode ")){
            	
            	Vector<String> gpv = new Vector<String>();
            	
                Pattern genename = Pattern.compile("TextLabel=\"(.*?)\".*Type=\"GeneProduct\"");//

                Matcher mg = genename.matcher(ln);
                
                String genam = "null";
                
                if(mg.find()){
                	
                	String genamP = mg.group(1);
                	
                	genam = genamP.replaceAll("\\?", "").replaceAll("`","").replaceAll("@", "");//.toUpperCase()                	
                	
                    gv.set(0,genam);

                    gpv.add(genam);

                    StoreGenes(genam, pathName);

                }

                Pattern graphID = Pattern.compile("GraphId=\"(.*?)\"");
                
                Matcher mgra = graphID.matcher(ln);
                
                while(mgra.find()){
                	
                    k = mgra.group(1);

                }

                Pattern groupID = Pattern.compile("GroupRef=\"(.*?)\"");
                
                Matcher mgp = groupID.matcher(ln);

                if(mgp.find()){
                	
                    gpk = mgp.group(1);

                    //System.out.println(gpk);
                    if(!groupmp.containsKey(gpk)){
                    	
                        groupmp.put(gpk,gpv);
                        
                        gpv = new Vector<String>();
                        
                    }else{
                        //if(mg.find()){
                    	
                    	if(genam!="null"){

                            //System.out.println(mg.group(1));
                            Vector<String> tmpv = groupmp.get(gpk);
                            
                            if(!tmpv.contains(genam)){
                            	
                            	tmpv.add(genam);
                            	
                            	genam = "null";
                            	
                            	groupmp.put(gpk,tmpv);
                            	
                            }
                                                                                    
                            //gpv = new Vector<String>();

                    	}
                    }
                    
                    gpv = new Vector<String>();
                    
                }
            }

            if(ln.contains("<Xref Database=")){
                
            	//This was a bug!
            	//Pattern geneid = Pattern.compile("<Xref Database=\"(.*)\" ID=\"(\\w*)\" />");
            	Pattern geneid = Pattern.compile("<Xref Database=\"(.*)\" ID=\"(.*)\"");
            	
                Matcher mgid = geneid.matcher(ln);
                
                while(mgid.find()){
                    //gv.set(1,mgid.group(1)+"\t"+mgid.group(2));
                	
                	gv.set(1,mgid.group(2));
                	
                }
            }
            
            if(ln.contains("</DataNode>")){
            	
            	if(k==null || gv.get(0)==null) continue;
                
            	if(gragenmp.containsKey(k)){
                    
            		System.out.println("Warning! Something Wrong here, in GPML two entry have same graph ID.");

                }else {
                	                	
                	gragenmp.put(k, gv);
                	
                }   
            	
        		//This wikimap is to used GPML inherent ID information to do further mapping.
        		
        		if(gv.get(1)!=null && gv.get(1).length()>1){
        			
        			wikimp.put(gv.get(0),gv.get(1));
        			
        			//System.out.println(gv); 
        		}
        		
            	
            	
                gv = new ArrayList<String>();
                
                gv.add(null); gv.add(null);                
            }
     
           
            

//            if(ln.contains("</DataNode>")){
//                
//            	if(k!=null && (gv.get(0))!=null && !gragenmp.containsKey(k)){
//                    
//                	gragenmp.put(k, gv);
//
//                }                
//                
//                gv = new ArrayList<String>();
//                
//                gv.add(null); gv.add(null);                
//            }

            if(ln.contains("<Point")){
                //Pattern nodp = Pattern.compile("GraphRef=\"(\\w*)\"");
                Pattern nodp = Pattern.compile("GraphRef=\"(.*?)\"");
                
                Matcher nodm = nodp.matcher(ln);
                
                while(nodm.find()){
                    //System.out.println(nodm.group(1));
                    
                	nodv.add(nodm.group(1));
                    //System.out.println("debug");
                	
                }
            }

            if(ln.contains("</Line>")){
            	
                if(nodv.size()==2){
                	
                	linepairv.add(nodv);
                	
                }
                
                nodv = new ArrayList<String>();
                
            }

        }
        br.close();
        //System.out.println(pathName+"\t"+groupmp.size());
    }
    
    private static void writepair(Vector<ArrayList<String>> linepairv, 
            HashMap<String, ArrayList> gragenmp, HashMap<String, Vector> groupmp, 
            String odir) throws IOException {

        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(odir)));
        
        HashMap<String,Boolean> ckamp = new HashMap<String,Boolean>();
        
        HashMap<String,Boolean> ckbmp = new HashMap<String,Boolean>();

        for(int i = 0; i<linepairv.size();i++){
            
        	ArrayList<String> eleal = linepairv.get(i);
            
        	String na = eleal.get(0);
            
        	String nb = eleal.get(1);
            
        	if(gragenmp.containsKey(na)&&gragenmp.containsKey(nb)){
                
        		ArrayList<String> genA = gragenmp.get(na);
                
        		ArrayList<String> genB = gragenmp.get(nb);
        		
//        		if(genA.get(1)!=null && genA.get(1).length()>1){
//        			
//        			wikimp.put(genA.get(0),genA.get(1));
//        		}
//        		
//        		if(genB.get(1)!=null && genB.get(1).length()>1){
//        			
//        			wikimp.put(genB.get(0),genB.get(1));
//        		}
        		
                //System.out.println(genA.get(0)+"\t"+"|"+"\t"+genB.get(0) +"\t"+"|"+"\t"+"graphRel"+"\t"+"|"+"\t"+pathName+"\t"+"|"+"\t"+genA.get(1)+"\t"+"|"+"\t"+genB.get(1));
                
        		if(genA.get(0).length()>0 && genB.get(0).length()>0){//genA.get(0).length()<10 && genB.get(0).length()<10 && 
                    
        			if(!ckamp.containsKey(genA.get(0)+"\t"+genB.get(0))&&!ckamp.containsKey(genB.get(0)+"\t"+genA.get(0))){
                                				                       
        				ckamp.put(genA.get(0)+"\t"+genB.get(0), Boolean.TRUE);
        				
        				String gepreA = genA.get(0);
        				
        				String geA = gepreA.replaceAll("\\?", "").replaceAll("`","").replaceAll("@", "");
        				
        				String gepreB = genB.get(0);
        				
        				String geB = gepreB.replaceAll("\\?", "").replaceAll("`","").replaceAll("@", "");        				
        				        				
        				pw.println(geA+"\t"+geB +"\t"+"graphRel"+"\t"+pathName+"\t");
                        
        				StoreGenePairs(geA, geB, "graphRel", pathName);
        				
                        //System.out.println("hello>"+genA.get(0)+"\t"+genB.get(0));
                    }
                }
            }
        }
        
        HashMap<String,Vector<String>> egromp = new HashMap<String,Vector<String>>();
        
        Set<String> groupsets = groupmp.keySet();        
        
        for(String id:groupsets){            
            
        	Vector<String> egrovv = new Vector<String>();
        	
        	Vector env = (Vector) groupmp.get(id);
            
        	for(int i = 0; i<env.size();i++){
                
    			String genea = (String) env.get(i);
                String gena = genea.replaceAll("\\?", "").replaceAll("`","");//replaceAll(" ", "").toUpperCase()
    			
                if(!egrovv.contains(gena)) {
                	
                	egrovv.add(gena);
                	
                	//System.out.println("debug: group container");
                }
                
                
        		
        		for(int j = i+1; j<env.size();j++){ 
                    
                    String geneb = (String) env.get(j); 
        			String genb = geneb.replaceAll("\\?", "").replaceAll("`","");//replaceAll(" ", "").toUpperCase()                            		
        			
        			//System.out.println(gena+"\t"+"|"+"\t"+genb +"\t"+"|"+"\t"+"groupRel"+"\t"+"|"+"\t"+pathName);
                    
        			if(gena.length()<10&&genb.length()<10&&!gena.equalsIgnoreCase(genb)){
                        
        				if(!ckbmp.containsKey(gena+"\t"+genb)&&!ckbmp.containsKey(genb+"\t"+gena)){
                            
        					ckbmp.put(gena+"\t"+genb, Boolean.TRUE);
                            
        					pw.println(gena+"\t"+genb +"\t"+"groupRel"+"\t"+pathName);
                            
        					StoreGenePairs(gena, genb, "groupRel", pathName);
        					
                        }                       
                    }
                }
            }
            
        	egromp.put(id, egrovv);
        	
        	//System.out.println("Size of Element Group HashMap: "+egromp.size());
        	
        }
        
        phgro.put(pathName, egromp);
        
        //System.out.println("Size of Group container: "+phgro.size());
        
        
//        Set<String> groupsets = groupmp.keySet();        
//        
//        for(String id:groupsets){            
//            
//        	Vector env = (Vector) groupmp.get(id);
//            
//        	for(int i = 0; i<env.size();i++){
//                
//        		for(int j = i+1; j<env.size();j++){ 
//                    
//        			String genea = (String) env.get(i);
//                    String gena = genea.replaceAll("\\?", "").replaceAll("`","");//replaceAll(" ", "").toUpperCase()
//        			
//                    String geneb = (String) env.get(j); 
//        			String genb = geneb.replaceAll("\\?", "").replaceAll("`","");//replaceAll(" ", "").toUpperCase()
//                    
//        			//System.out.println(gena+"\t"+"|"+"\t"+genb +"\t"+"|"+"\t"+"groupRel"+"\t"+"|"+"\t"+pathName);
//                    
//        			if(gena.length()<10&&genb.length()<10&&!gena.equalsIgnoreCase(genb)){
//                        
//        				if(!ckbmp.containsKey(gena+"\t"+genb)&&!ckbmp.containsKey(genb+"\t"+gena)){
//                            
//        					ckbmp.put(gena+"\t"+genb, Boolean.TRUE);
//                            
//        					pw.println(gena+"\t"+genb +"\t"+"groupRel"+"\t"+pathName);
//                            
//        					StoreGenePairs(gena, genb, "groupRel", pathName);
//        					
//                        }                       
//                    }
//                }
//            }
//            
//        }
        
        pw.close();
        
    }	
    
    private static void StoreGenes(String geneName,String pathway) {
        
    	HashMap<String, String> ges = new HashMap<String,String>();  
        
        //if (arr[0].length()<1) continue;
        
        if (phge.containsKey(pathway)){
        	
        	ges = phge.get(pathway);
            
        	if (!ges.containsKey(geneName)){
            	
        		ges.put(geneName,"WikiPathWays");
        		
        		phge.put(pathway, ges);        		
            	
            }            	        	        	
        	
        } else {
        	
        	ges.put(geneName, "WikiPathways");
        	
        	phge.put(pathway, ges);
        	
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
    
    public static void main(String args[]) throws IOException{
    
    	extract("sapiens");
    	
    }
	
}
