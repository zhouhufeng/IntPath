/* Title: Biocyc parser
 * Function: Java Program to Extrace Information from OWL 
 * 
 * Author: Hufeng Zhou
 * Time: July 3rd 2021
 * Version: v2
 */


package sapiens;

import java.io.*;
import org.biopax.paxtools.controller.EditorMap;
import org.biopax.paxtools.io.BioPAXIOHandler;
import org.biopax.paxtools.io.sif.SimpleInteractionConverter;
import org.biopax.paxtools.io.sif.level2.ComponentRule;
import org.biopax.paxtools.io.sif.level2.ConsecutiveCatalysisRule;
import org.biopax.paxtools.io.sif.level2.ControlRule;
import org.biopax.paxtools.io.sif.level2.ControlsTogetherRule;
import org.biopax.paxtools.io.sif.level2.ParticipatesRule;
import org.biopax.paxtools.io.simpleIO.SimpleEditorMap;
import org.biopax.paxtools.io.simpleIO.SimpleReader;
import org.biopax.paxtools.model.*;
//import org.biopax.paxtools.io.sif.level3.ComponentRule;
//import org.biopax.paxtools.io.sif.level3.ConsecutiveCatalysisRule;
//import org.biopax.paxtools.io.sif.level3.ControlRule;
//import org.biopax.paxtools.io.sif.level3.ControlsTogetherRule;
//import org.biopax.paxtools.io.sif.level3.ParticipatesRule;
import java.util.*;
import java.util.regex.*;
import tools.*;

public class BioCyc {
	
	//The HashMap is used to store the information of which genes are in certain pathways. With pathway name as the key, gene names stored in a HashMap as a value. 
	public static HashMap<String,HashMap<String,String>> phge = new HashMap<String,HashMap<String,String>>(10000);
	/* The HashMap is used to store the information of which gene pairs are in certain pathways. With pathway 
	 * name as the key, Value is a HashMap, gene pair names (A>B) as key stored , relationship as value in 
	 * this HashMap as a value.
	 */
	public static HashMap<String,HashMap<String,String>> phgp = new HashMap<String,HashMap<String,String>>(10000);
	
    public static void extract(String orgs) throws Exception {
    //public static void main(String args[]) throws IOException {
    	
		/*
		 * Convert the biopax file into sif file.
		 */
        String ibx = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"biopax-level2.owl";
        String nod = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"Humancycl2node.sif";
        String edg = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"Humancycl2edge.sif";
     	
        convertBiopaxToSIF(ibx,nod,edg);
   
            
		/*
		 * Maping the sif file nodes and edges, which is mapping the results of paxtools to get the specific result I want.
		 */
        String oclean = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"Humancyc2clearned.txt";
        String pair   = orgs+File.separator+"extraction"+File.separator+"BioCyc"+File.separator+"Humancyc2pair.txt";
        String nmp    = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"protein-links.dat";
                         
        HashMap<String,Vector<String>> ens = frmENS(nmp);
        
        HashMap nodmp = constructProNamMapping(nod,ens);
        
        writePairWiseResult(nodmp,edg,oclean,pair,ens);        
        
		/*
		 * The following step is used to add the pathway names to each of the extracted gene pairs.        
		 */
        String gepath = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"pathways.col";
        
        String gpair  = orgs+File.separator+"extraction"+File.separator+"BioCyc"+File.separator+"HumanCycPathPairs";
        
        Normalize.genecolmp(orgs);
        
        HashMap<String,Vector<String>>idmp = Normalize.mp;
        
        //test tidmp = new test();        
        
        //tidmp.writeOutVectorKey(idmp);
        
        HashMap<String, Vector<String>> pathmp = Constructpgmp(gepath);
        
        //tidmp.writeOutVectorKey(pathmp);

        Writepairrelwithpathname(pathmp,pair,gpair,idmp);        
        
//        //Construct two HashMap, one to map the MGI to the official symbol, 
//        HashMap<String,String> idmp = new HashMap<String,String>();
//        //the other HashMap to map the information of MGI as Key, and all the involving pathway name stored in a vector use as retrieval.
//        HashMap<String,Vector> pathmp = new HashMap<String,Vector>();
//                       
//        //This step is used to construct pathway genes map. 
//        Constructmp(idmp,pathmp,gepath);
//        
//        Writepairrelwithpathname(idmp,pathmp,pair,gpair);        
     }
	
    private static void convertBiopaxToSIF(String ibx, String nod, String edg)throws IOException {
		// TODO Auto-generated method stub
                BioPAXIOHandler biopaxReader = new SimpleReader();
                Model model = biopaxReader.convertFromOWL(new FileInputStream(new File(ibx)));
                SimpleInteractionConverter sic = null;
                if (BioPAXLevel.L2.equals(model.getLevel())) {
                        sic = new SimpleInteractionConverter(new ComponentRule(),
                                        new ConsecutiveCatalysisRule(), new ControlRule(),
                                        new ControlsTogetherRule(), new ParticipatesRule());
                } else if (BioPAXLevel.L2.equals(model.getLevel())) {
                        sic = new SimpleInteractionConverter(
                                        new org.biopax.paxtools.io.sif.level3.ComponentRule(),
                                        new org.biopax.paxtools.io.sif.level3.ConsecutiveCatalysisRule(),
                                        new org.biopax.paxtools.io.sif.level3.ControlRule(),
                                        new org.biopax.paxtools.io.sif.level3.ControlsTogetherRule(),
                                        new org.biopax.paxtools.io.sif.level3.ParticipatesRule());
                } else {
                        System.err.println("SIF converter does not yet support BioPAX level: "
                                        + model.getLevel());
                        System.exit(0);
                }

                EditorMap editorMap = new SimpleEditorMap(model.getLevel());
                OutputStream edgeStream = new FileOutputStream(new File(edg));
                OutputStream nodeStream = new FileOutputStream(new File(nod));
                sic.writeInteractionsInSIFNX(model, edgeStream, nodeStream,
                        false, editorMap,"NAME","XREF" ); //Bug identified,
            //for BioPAX L3 use, "name","xref"
            //for Biopax L2 use, "NAME","XREF"        		
	}

    
    private static HashMap<String, Vector<String>> frmENS (String nmp) throws IOException{
        
        String ln;
        
        String[] arr;
        
        HashMap<String,Vector<String>> idmp = new HashMap<String,Vector<String>>();
        
        BufferedReader br = new BufferedReader(new FileReader(nmp));
        
        while((ln = br.readLine())!=null){

            arr = ln.split("\t");
            
            if (arr.length < 2) continue;
            
            if(!ln.contains("ENS")) continue;
            
            String[] ens = arr[0].split("-");

            Vector<String> val = new Vector<String>();

            if(!idmp.containsKey(ens[0])){
            	
                val.add(arr[1]);
                
                idmp.put(ens[0],val);
                
            }else{
                
            	val = idmp.get(ens[0]);
                
            	if(!val.contains(arr[1])){
                    
            		val.add(arr[1]);
                    
            		idmp.put(ens[0],val);
            		
                }
            }            
        }
        
        return idmp;
    }      
        
    private static HashMap constructProNamMapping(String inod,HashMap<String,Vector<String>> ens) throws IOException{
        HashMap<String,String> nodmp = new HashMap<String,String>();
        
        BufferedReader br = new BufferedReader(new FileReader(inod));
        
        String ln ;
        
        String[] arr;

        while((ln = br.readLine())!=null){
            arr = ln.split("\t");
            String[] proar = arr[0].split("#");

            if(arr.length>2 && arr[2].contains("HumanCyc:")){
                String gene = "";
                Pattern mp = Pattern.compile("HumanCyc:(.*?)-MONOMER;");
                Matcher mm = mp.matcher(arr[2]);
                if(mm.find()){
                    //System.out.print(mm.group(1)+"\n");
                    //gene = mm.group(1);
                    String tmp = mm.group(1);
                    if(ens.containsKey(tmp)){
                        Vector<String> tv = ens.get(tmp);
                        for(int i = 0; i<tv.size();i++){
                            gene = tv.get(i);
                            nodmp.put(proar[1],gene);
                            //System.out.println("debug+"+gene);
                        }
                    }else{
                        gene = tmp;
                        nodmp.put(proar[1],gene);
                    }
                    //System.out.print(gene+"\n");
                }
//                if(gene.length() < 1) continue;
//                nodmp.put(proar[1],gene);
            }            
        }
        
        //System.out.println(nodmp.size());
        
        return nodmp;
    }
	    
    private static void writePairWiseResult(HashMap nodmp, String iedg, String ocleaned, String omapped, HashMap<String,Vector<String>> ens) throws IOException {
    
    	
        BufferedReader bredg = new BufferedReader(new FileReader(iedg));
        
        PrintWriter pwcl     = new PrintWriter(new BufferedWriter(new FileWriter(ocleaned)));
        
        PrintWriter pw 		 = new PrintWriter(new BufferedWriter(new FileWriter(omapped)));
        
        String ln;
        
        String[] arr;

        HashMap<String,Boolean> ckmp = new HashMap<String,Boolean>();

        while((ln = bredg.readLine())!=null){
        	
            arr = ln.split("\t");
            
            String[] ita = arr[0].split("#");
            
            String[] itb = arr[2].split("#");
            
            if(!ckmp.containsKey(ita[1]+"\t"+itb[1])&&!ckmp.containsKey(itb[1]+"\t"+ita[1])){
                
            	pwcl.println(ita[1]+"\t"+itb[1]+"\t"+arr[1]);
                
            	//System.out.println(ita[1]+"\t"+itb[1]+"\t"+arr[1]);
                if(nodmp.containsKey(ita[1])&&nodmp.containsKey(itb[1])){
                    
                	String va = (String) nodmp.get(ita[1]);
                    
                	String vb = (String) nodmp.get(itb[1]);
                    
                	if (!nodmp.containsKey(ita[1])||!nodmp.containsKey(itb[1])) continue;
                    
                	if (va.length()<1||vb.length()<1) continue;
                    
                	pw.println(va+"\t"+vb+"\t"+arr[1]+"\t");
                	
                }
                
                ckmp.put(ita[1]+"\t"+itb[1],Boolean.TRUE);
            }

        }
        
        bredg.close();
        
        pwcl.close();
        
        pw.close();
    }    
    
    private static HashMap Constructpgmp(String pgm) throws IOException {

        HashMap<String, Vector> pathmp = new HashMap<String, Vector>();
        
        BufferedReader br = new BufferedReader(new FileReader(pgm));

        String ln, id;
        
        String[] arr;

        while ((ln = br.readLine()) != null) {
            
        	arr = ln.split("\t");
            
            for (int i = 2; i < arr.length && i<95; i++) {
                //if(pidmp.containsKey(arr[0])){
            		
            		if(arr[i].length()<2) break;
                    
            		Vector<String> pv = new Vector<String>();
                    
                    if (!pathmp.containsKey(arr[i])) {
                    	
                        pv.add(arr[1]);
                        
                        pathmp.put(arr[i], pv);                                                
                        
                    } else {
                        
                    	pv = pathmp.get(arr[i]);
                        
                        if (!pv.contains(arr[1])) {
                        
                        	pv.add(arr[1]);
                        	
                        }
                        
                        pathmp.put(arr[i], pv);
                    }
                    
                  StoreGenes(arr[i],arr[1]);
                  
                  HashMap<String,String> genes = new HashMap<String,String>();
                  
//                  if(phge.containsKey(arr[1])){
//                  	
//                	genes = phge.get(arr[1]);
//                	
//                  	if (!genes.containsKey(arr[i])){            		            		
//                  		genes.put(arr[i],"BioCyc");
//                      	phge.put(arr[1],genes);
//                  	}
//                  	
//                  }else{
//                	  
//                  	genes.put(arr[i],"BioCyc");
//                  	phge.put(arr[1],genes);
//                  	
//                  }                    
                    
                //}
            }
        }
        br.close();
        
        //System.out.println("The size of pathway-genes map is: " + pathmp.size());
        
        return pathmp;
    }

    private static void Writepairrelwithpathname(HashMap<String, Vector<String>> pathmp, String pair, String res,
    		HashMap<String,Vector<String>>idmp) throws IOException {

        String ln;
        
        String[] arr;

        BufferedReader br = new BufferedReader(new FileReader(pair));
        
        PrintWriter    pw = new PrintWriter(new BufferedWriter(new FileWriter(res)));

        HashMap<String, Boolean> ckmp = new HashMap<String, Boolean>();

        while ((ln = br.readLine()) != null) {
        	
            arr = ln.split("\t");
            
            if(arr[0].length()<1||arr[1].length()<1) continue;
            
            Vector<String> gva = new Vector<String>();
            
            Vector<String> gvb = new Vector<String>();
            
            if(idmp.containsKey(arr[0])&&idmp.containsKey(arr[1])){
            	
            	gva = idmp.get(arr[0]);
            	
            	gvb = idmp.get(arr[1]);
            	
            }else if(idmp.containsKey(arr[0])&&!idmp.containsKey(arr[1])){
            	
            	gva = idmp.get(arr[0]);
            	
            	gvb.add(arr[1]);                                    	
            	
            }else if(!idmp.containsKey(arr[0])&&idmp.containsKey(arr[1])){
            	
            	gva.add(arr[0]);
            	
            	gvb = idmp.get(arr[1]);                                  	
            	
            }else if(!idmp.containsKey(arr[0])&&idmp.containsKey(arr[1])){
            	
            	gva.add(arr[0]);
            	
            	gvb.add(arr[1]);                                   	
            	
            }
            
            
            
            for(int p = 0; p< gva.size();p++){
            	
            	for(int q = 0; q<gvb.size();q++){
            		
            		if(gva.get(p).length()<2||gvb.get(q).length()<2) continue;
            		
            		String gp = gva.get(p);
            		
            		String gq = gvb.get(q);
            		            		            
//            		String gp = gva.get(p).replaceAll("\\?", "").replaceAll(" ", "").replaceAll("`","");
//            		
//            		String gq = gvb.get(q).replaceAll("\\?", "").replaceAll(" ", "").replaceAll("`","");
//            		

                    
                    if (pathmp.containsKey(gp) && pathmp.containsKey(gq)) {

                        //System.out.println("debug"+arr[0]+arr[0].toUpperCase()+arr[1]+arr[1].toUpperCase());
                        Vector pva = pathmp.get(gp);
                        
                        Vector pvb = pathmp.get(gq);
                        //System.out.println(pva.size()+"\t"+pvb.size());
                        Vector<String> ckv = new Vector<String>();
                        //if(pva.size()<2||pvb.size()<2) continue;
                        
                        for (int a = 0; a < pva.size(); a++) {
                            
                        	for (int b = 0; b < pvb.size(); b++) {
                            	
                                if (pva.get(a).equals(pvb.get(b))) {  
                                	
                            		if(!ckmp.containsKey(gp + "\t" + gq + "\t" + pva.get(a))){    
                        			
                                    pw.println(gp + "\t" + gq + "\t" + arr[2] + "\t" + pva.get(a));
                                    
                                    ckmp.put(gp + "\t" + gq + "\t" + pva.get(a), Boolean.TRUE);
                                    
                        		}
                                
                                StoreGenePairs(gp, gq, arr[2], (String)pva.get(a));                                	
                                	
                                            
//                            		if(!ckmp.containsKey(arr[0] + "\t" + arr[1] + "\t" + pva.get(a))){    
//                            			
//                                        pw.println(arr[0] + "\t" + arr[1] + "\t" + arr[2] + "\t" + pva.get(a));
//                                        
//                                        ckmp.put(arr[0] + "\t" + arr[1] + "\t" + pva.get(a), Boolean.TRUE);
//                                        
//                            		}
//                                    
//                                    StoreGenePairs(arr[0], arr[1], arr[2], (String)pva.get(a));                       
                                	
                                }                                               
                            }
                        }
                    }                    
                    
                    
            		
            	}
            	
            }            
            
            
            

        }
        br.close();
        pw.close();
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
    
//    private static void Constructmp(HashMap<String, String> idmp, HashMap<String, Vector> pathmp, String gepath) throws IOException{
//    	
//        BufferedReader br = new BufferedReader(new FileReader(gepath));
//        
//        String ln, id;
//        
//        String[] arr;
//        
//        while((ln = br.readLine())!=null){
//            
//            Vector<String> pv = new Vector<String>();
//            
//            arr = ln.split("\t");
//            
//            if (arr.length<4) continue;
//            
//            HashMap<String,String> genes = new HashMap<String,String>();
//            
//            if(phge.containsKey(arr[3])){
//            	genes = phge.get(arr[3]);
//            	if (!genes.containsKey(arr[1])){            		            		
//            		genes.put(arr[1],"BioCyc");
//            	}            
//            	phge.put(arr[3],genes);
//            	
//            }else{
//            	
//            	genes.put(arr[1],"BioCyc");
//            	phge.put(arr[3],genes);
//            }
//            
//            Pattern idp = Pattern.compile("MGI:(\\d{3,15})");
//            Matcher idm = idp.matcher(ln);
//            
//            if(idm.find()){
//                id = idm.group(1);
//                //System.out.println(id);            
//            
//                if(!idmp.containsKey(id)){
//                    idmp.put(id, arr[1]);
//                }
//                arr[3] = arr[3].trim();
//                if(!pathmp.containsKey(id)){  
//                    pv.add(arr[3]);
//                    pathmp.put(id,pv);                
//                }else{
//                    pv = pathmp.get(id);
//                    if(!pv.contains(arr[3])){
//                        pv.add(arr[3]);
//                    }
//                    pathmp.put(id,pv);
//                }
//            }
//        }
//        br.close();
//        System.out.println( "The size of node ID map is: "+idmp.size()+"\n"+ "The size of node pathway name map is: "+pathmp.size());
//        
//    }
//
//    private static void Writepairrelwithpathname(HashMap<String, String> idmp, HashMap<String, Vector> pathmp, String pair, String res)throws IOException {
//        
//        String ln;
//        String[] arr;        
//        
//        BufferedReader br = new BufferedReader(new FileReader(pair));
//        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(res)));
//        
//        HashMap<String,Boolean> ckmp = new HashMap<String,Boolean>();
//
//        while((ln = br.readLine())!=null){
//            
//        	arr = ln.split("\t");
//            
//            
//            if(pathmp.containsKey(arr[6]) && pathmp.containsKey(arr[7])){
//                
//                //System.out.println("debug");
//                Vector pva = pathmp.get(arr[6]);
//                
//                Vector pvb = pathmp.get(arr[7]);
//                
//                Vector<String> ckv = new Vector<String>();
//                
//                for(int a = 0;a<pva.size();a++){
//                    
//                	for(int b = 0;b <pvb.size();b++){
//                    
//                		if(pva.get(a).equals(pvb.get(b))){
//                        
//                			if(!ckv.contains(((String) pva.get(a)).trim())){
//                                
//                            	ckv.add(((String) pva.get(a)).trim());
//                                
//                                String geneA = idmp.get(arr[6]);
//                                
//                                String geneB = idmp.get(arr[7]);
//                                //if(!ckmp.containsKey(geneA+"\t"+geneB)&&!ckmp.containsKey(geneB+"\t"+geneA)){
//                                	
//                                	String pathway = (String) pva.get(a);
//                                    
//                                	pw.println(geneA+"\t"+geneB+"\t"+arr[5]+"\t"+pathway);
//                                    
//                                	ckmp.put(geneA+"\t"+geneB, Boolean.TRUE);                                                                        
//                                    
//                                    HashMap<String,String> pairs = new HashMap<String,String>();
//                                    
//                                    String ga, gb;
//                                    
//                                    if(geneA.compareTo(geneB)>=0){
//                                    	
//                                    	ga = geneA;
//                                    	
//                                    	gb = geneB;
//                                    	
//                                    }else{
//                                    	
//                                    	gb = geneA;
//                                    	
//                                    	ga = geneB;
//                                    	
//                                    }
//                                    
//                                    if(phgp.containsKey(pva.get(a))){
//                                    	
//                                    	pairs = phgp.get(pathway);
//                                    	
//                                    	if(!pairs.containsKey(ga+"\t"+gb)){
//                                    		
//                                    		pairs.put(ga+"\t"+gb, arr[5]);
//                                    	}
//                                    	
//                                    	phgp.put(pathway, pairs);
//                                    	
//                                    }else{
//                                    	
//                                    	pairs.put(ga+"\t"+gb, arr[5]);
//                                    	
//                                    	phgp.put(pathway, pairs);
//                                    	
//                                    }
//                                //}
//                            }                           
//                        }
//                    }
//                }                
//            }
//        }
//        
//        br.close();
//        
//        pw.close();
//    }

    
    private static void StoreGenes(String geneName,String pathway) {
        
    	HashMap<String, String> ges = new HashMap<String,String>();  
        
        //if (arr[0].length()<1) continue;
        
        if (phge.containsKey(pathway)){
        	
        	ges = phge.get(pathway);
            
        	if (!ges.containsKey(geneName)){
            	
        		ges.put(geneName,"BioCyc");
        		
        		phge.put(pathway, ges);        		
            	
            }            	        	        	
        	
        } else {
        	
        	ges.put(geneName, "BioCyc");
        	
        	phge.put(pathway, ges);
        	
        }  
		
	}    
}
