package tuberculosis;

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
	
	
	public static void main(String args[])throws Exception{
		
		extract("tuberculosis");
		
	}
	
	
    public static void extract(String orgs) throws Exception {
    //public static void main(String args[]) throws IOException {
    	
		/*
		 * Convert the biopax file into sif file.
		 */
        String ibx = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"biopax-level2.owl";
        String nod = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+orgs+"Cycl2node.sif";
        String edg = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+orgs+"Cycl2edge.sif";
     	
        convertBiopaxToSIF(ibx,nod,edg);
               
		/*
		 * Maping the sif file nodes and edges, which is mapping the results of paxtools to get the specific result I want.
		 */
        String oclean = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+orgs+"Cyc2clearned.txt";
        
        String pair   = orgs+File.separator+"extraction"+File.separator+"BioCyc"+File.separator+orgs+"Cyc2pair.txt";
        
        HashMap<String, Vector<String>> nodmp = constructProNamMapping(nod);
        
        writePairWiseResult(nodmp,edg,oclean,pair);        
        
		/*
		 * The following step is used to add the pathway names to each of the extracted gene pairs.        
		 */
        String gepath = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"pathways.col";
        
        String gpair  = orgs+File.separator+"extraction"+File.separator+"BioCyc"+File.separator+orgs+"CycPathPairs";
        
        //Normalize.genecolmp(orgs);
                
        HashMap<String, Vector<String>> pathmp = Constructpgmp(gepath);
        
        Writepairrelwithpathname(pathmp,pair,gpair);        
           
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

                EditorMap editorMap 	= new SimpleEditorMap(model.getLevel());
                OutputStream edgeStream = new FileOutputStream(new File(edg));
                OutputStream nodeStream = new FileOutputStream(new File(nod));
                sic.writeInteractionsInSIFNX(model, edgeStream, nodeStream,
                        false, editorMap,"NAME","XREF" ); //Bug identified,       		
	}
   
    private static HashMap constructProNamMapping(String inod) throws IOException{
        
    	HashMap<String,Vector<String>> nodmp = new HashMap<String,Vector<String>>();
        
        BufferedReader br = new BufferedReader(new FileReader(inod));
        
        String ln ;
        
        String[] arr;

        while((ln = br.readLine())!=null){
            
        	arr = ln.split("\t");
            
        	String[] proar = arr[0].split("#");

            if(arr.length>2 && arr[2].contains("MTBRvCyc:")){
            	
                String gene = "";
                
                Pattern mp = Pattern.compile("MTBRvCyc:(.*?)-MONOMER;");
                
                Matcher mm = mp.matcher(arr[2]);                
                
                if(mm.find()){

                    String ygen = mm.group(1);
                    
                    if(nodmp.containsKey(proar[1])){
                    	
                    	Vector<String> vv = nodmp.get(proar[1]);
                    	
                    	if(!vv.contains(ygen)){
                    		
                    		vv.add(ygen);
                    		
                    		nodmp.put(proar[1], vv);
                    		
                    	}
                    	
                    }else{
                    	
                    	Vector<String> vv = new Vector<String>();
                    	
                    	vv.add(ygen);
                    	
                    	nodmp.put(proar[1], vv);
                    }                    
                }
            }            
        }
        return nodmp;
    }
	    
    private static void writePairWiseResult(HashMap<String,Vector<String>> nodmp, String iedg, 
    		String ocleaned, String omapped) throws IOException {
    
    	
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
                	
                	//System.out.println("Debug! BioCyc extraction, map pair to pathways");
                	
                	Vector<String> vva = nodmp.get(ita[1]);
                	
                	Vector<String> vvb = nodmp.get(itb[1]);
                	
                	for(int a = 0; a< vva.size();a++){
                		
                		for(int b = 0; b<vvb.size();b++){
                	
                        	String va = vva.get(a);
                            
                        	String vb = vvb.get(b);
                            
                        	if (!nodmp.containsKey(ita[1])||!nodmp.containsKey(itb[1])) continue;
                            
                        	if (va.length()<1||vb.length()<1) continue;
                            
                        	pw.println(va+"\t"+vb+"\t"+arr[1]+"\t");
                			
                		}                		
                		
                	}
                    
                	
                }
                
                ckmp.put(ita[1]+"\t"+itb[1],Boolean.TRUE);
            }

        }
        
        bredg.close();
        
        pwcl.close();
        
        pw.close();
    }    
    
    private static HashMap Constructpgmp(String pgm) throws IOException {

        HashMap<String, Vector<String>> pathmp = new HashMap<String, Vector<String>>();
        
        BufferedReader br = new BufferedReader(new FileReader(pgm));

        String ln, id;
        
        String[] arr;

        while ((ln = br.readLine()) != null) {
            
        	arr = ln.split("\t");
            
            for (int i = 48; i < arr.length && i<95; i++) {
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
                  
            }
        }
        br.close();
        
        //System.out.println("The size of pathway-genes map is: " + pathmp.size());
        
        return pathmp;
    }
   
    

    
    private static void Writepairrelwithpathname(HashMap<String, Vector<String>> pathmp, 
    		String pair, String res) throws IOException {

        String ln;
        
        String[] arr;

        BufferedReader br = new BufferedReader(new FileReader(pair));
        
        PrintWriter    pw = new PrintWriter(new BufferedWriter(new FileWriter(res)));

        HashMap<String, Boolean> ckmp = new HashMap<String, Boolean>();

        while ((ln = br.readLine()) != null) {
        	
            arr = ln.split("\t");
            
            if(arr[0].length()<1||arr[1].length()<1) continue;
                       		
    		String gp = arr[0];
    		
    		String gq = arr[1];            		            		                    		

            
            if (pathmp.containsKey(gp) && pathmp.containsKey(gq)) {

                //System.out.println("debug"+arr[0]+arr[0].toUpperCase()+arr[1]+arr[1].toUpperCase());
                Vector pva = pathmp.get(gp);
                
                Vector pvb = pathmp.get(gq);
                
                Vector<String> ckv = new Vector<String>();
                
                for (int a = 0; a < pva.size(); a++) {
                    
                	for (int b = 0; b < pvb.size(); b++) {
                    	
                        if (pva.get(a).equals(pvb.get(b))) {  
                        	
                    		if(!ckmp.containsKey(gp + "\t" + gq + "\t" + pva.get(a))){    
                			
                            pw.println(RvIDcorrect(gp) + "\t" + RvIDcorrect(gq) + "\t" + arr[2] + "\t" + pva.get(a));
                            
                            ckmp.put(gp + "\t" + gq + "\t" + pva.get(a), Boolean.TRUE);
                            
                		}
                        
                        StoreGenePairs(gp, gq, arr[2], (String)pva.get(a));                                	                                	                                                             
                        	
                        }                                               
                    }
                }
            }                    
                    
        }
        
        br.close();
        
        pw.close();
    }
        
    private static void StoreGenePairs(String GeneA, String GeneB, String rel,String pathway) {
		
    	String geneA = RvIDcorrect(GeneA);
    	
    	String geneB = RvIDcorrect(GeneB);
    	
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
    

    
    private static String RvIDcorrect(String ins) {
		
    	String ots = "";
		
		if(ins.contains("RV")){
			
            String tmp = ins;
            String tmpe = "R"+tmp.substring(1,tmp.length()).toLowerCase();
            if(tmp.length()==6){
            	ots = tmpe;
            }else if(tmp.length()>6){
                if(tmp.substring(6,tmp.length()).equalsIgnoreCase("c")){
                	ots = tmpe;
                }else{
                    //System.out.print("bugs bugs\n");
                	ots = "R"+tmp.substring(1,6).toLowerCase()+tmp.substring(6,tmp.length());
                }
            }

			
		}else{
			
			ots = ins;
			
		}
		
		return ots;
	}

	private static void StoreGenes(String GeneName,String pathway) {
        
		String geneName = RvIDcorrect(GeneName);
		
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
