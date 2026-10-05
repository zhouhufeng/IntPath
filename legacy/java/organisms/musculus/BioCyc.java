package musculus;

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

public class BioCyc {
	
	//The HashMap is used to store the information of which genes are in certain pathways. With pathway name as the key, gene names stored in a HashMap as a value. 
	public static HashMap<String,HashMap<String,String>> phge = new HashMap<String,HashMap<String,String>>();
	/* The HashMap is used to store the information of which gene pairs are in certain pathways. With pathway 
	 * name as the key, Value is a HashMap, gene pair names (A>B) as key stored , relationship as value in 
	 * this HashMap as a value.
	 */
	public static HashMap<String,HashMap<String,String>> phgp = new HashMap<String,HashMap<String,String>>();
	
    public static void extract(String orgs) throws IOException {
    //public static void main(String args[]) throws IOException {
    	
		/*
		 * Convert the biopax file into sif file.
		 */
        String ibx = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"biopax-level2.owl";
        String nod = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"mousecycl2node.sif";
        String edg = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"mousecycl2edge.sif";
     	
        convertBiopaxToSIF(ibx,nod,edg);
   
            
		/*
		 * Maping the sif file nodes and edges, which is mapping the results of paxtools to get the specific result I want.
		 */
        String oclean = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"mousecyc2clearned.txt";
        String pair   = orgs+File.separator+"extraction"+File.separator+"BioCyc"+File.separator+"mousecyc2pair.txt";
        
        HashMap nodmp = constructProNamMapping(nod);
        writePairWiseResult(nodmp,edg,oclean,pair);             
        
        
		/*
		 * The following step is used to add the pathway names to each of the extracted gene pairs.        
		 */
        String gepath = orgs+File.separator+"source"+File.separator+"BioCyc"+File.separator+"pathways_summary.txt";
        
        String gpair  = orgs+File.separator+"extraction"+File.separator+"BioCyc"+File.separator+"MouseCycPairs";
        
        //Construct two HashMap, one to map the MGI to the official symbol, 
        HashMap<String,String> idmp = new HashMap<String,String>();
        //the other HashMap to map the information of MGI as Key, and all the involving pathway name stored in a vector use as retrieval.
        HashMap<String,Vector> pathmp = new HashMap<String,Vector>();                       
        //This step is used to
        Constructmp(idmp,pathmp,gepath);
        
        Writepairrelwithpathname(idmp,pathmp,pair,gpair);        
     }
	
    private static void convertBiopaxToSIF(String ibx, String nod, String edg)throws IOException {
		// TODO Auto-generated method stub
        
        // import BioPAX from OWL file (auto-detects level)
 		BioPAXIOHandler biopaxReader = new SimpleReader();
 		Model model = biopaxReader.convertFromOWL(new FileInputStream(new File(ibx)));
	    SimpleInteractionConverter sic = null;
 		if (BioPAXLevel.L2.equals(model.getLevel())) {
 			sic = new SimpleInteractionConverter(new ComponentRule(),
 					new ConsecutiveCatalysisRule(), new ControlRule(),
 					new ControlsTogetherRule(), new ParticipatesRule());
 		} else if (BioPAXLevel.L3.equals(model.getLevel())) {
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

	private static HashMap constructProNamMapping(String inod) throws IOException{
		
        HashMap<String,String[]> nodmp = new HashMap<String,String[]>();
        
        BufferedReader br = new BufferedReader(new FileReader(inod));
        
        String ln ;
        
        String[] arr;
        
        while((ln = br.readLine())!=null){
            arr = ln.split("\t");
            String[] proar = arr[0].split("#");
            String[] val = new String[4];
            val[0] = arr[1];
            val[3] = "NULL";
            //System.out.print(val[0]+"\t"+val[1]+"\t");
            if(arr.length>2 && arr[2].contains("UniProt:")){
                Pattern up = Pattern.compile("UniProt:(\\w*?);");
                Matcher um = up.matcher(arr[2]);
                if(um.find()){
                    //System.out.print(um.group(1)+"\t");
                    val[1]=um.group(1);
                }
            }else if(ln.contains("MGI:")){
                    val[1]="NULL";
                    //System.out.print(val[1]+"\t");
            }

            if(arr.length>2 && arr[2].contains("MGI:")){//MouseCyc:MGI:3026877-MONOMER;
                Pattern mp = Pattern.compile("MGI:(.*?)-MONOMER;");
                Matcher mm = mp.matcher(arr[2]);
                if(mm.find()){
                    //System.out.print(mm.group(1)+"\t");
                    val[2]=mm.group(1);
                }
            }   
            
            if(arr.length>2 && arr[2].contains("LIGAND-compound:")){
                Pattern up = Pattern.compile("LIGAND-compound:(\\w*?);");
                Matcher um = up.matcher(arr[2]);
                if(um.find()){
                    //System.out.print(um.group(1)+"\t");
                    val[1]=um.group(1);
                }
            }else if(ln.contains("MouseCyc:")&&!arr[2].contains("MGI:")){
                    val[1]="NULL";
                    //System.out.print(val[1]+"\t");
            }            

            if(arr.length>2 && arr[2].contains("MouseCyc:")&&!arr[2].contains("MGI:")){
                Pattern mp = Pattern.compile("MouseCyc:(.*?);");
                Matcher mm = mp.matcher(arr[2]);
                if(mm.find()){
                    //System.out.print(mm.group(1)+"\t");
                    //val[2]=mm.group(1);
                	val[2]="null";
                }                              
            }    
            
            if(arr.length>2 && ln.contains("smallMolecule")){
                Pattern cp = Pattern.compile("PubChem:(\\d*?);");
                Matcher cm = cp.matcher(arr[2]);
                if(cm.find()){
                    val[3]=cm.group(1);
                    //System.out.println(cm.group(1)+"\t");
                }
            }
            nodmp.put(proar[1],val);
            //System.out.println();
            
        }
        
        return nodmp;
    }

    private static void writePairWiseResult(HashMap nodmp, String iedg, String ocleaned, String omapped) throws IOException {
        BufferedReader bredg = new BufferedReader(new FileReader(iedg));
        PrintWriter pwcl = new PrintWriter(new BufferedWriter(new FileWriter(ocleaned)));
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(omapped)));
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
                    String[] va = (String[]) nodmp.get(ita[1]);
                    String[] vb = (String[]) nodmp.get(itb[1]);
                    pw.println(va[0]+"\t"+vb[0]+"\t"+arr[1]+"\t"+va[1]+"\t"+vb[1]+"\t"+arr[1]+"\t"+va[2]+"\t"+vb[2]);
                }                                         
                ckmp.put(ita[1]+"\t"+itb[1],Boolean.TRUE);
            }           
   
        }
        bredg.close();
        pwcl.close();
        pw.close();
    }	

    private static void Constructmp(HashMap<String, String> idmp, HashMap<String, Vector> pathmp, String gepath) throws IOException{
    	
        BufferedReader br = new BufferedReader(new FileReader(gepath));
        
        String ln, id;
        
        String[] arr;
        
        while((ln = br.readLine())!=null){
            
            
            
            arr = ln.split("\t");
            
            if (arr.length<4) continue;
            
            HashMap<String,String> genes = new HashMap<String,String>();
            
            if(phge.containsKey(arr[3])){
            	
            	genes = phge.get(arr[3]);
            	
            	if (!genes.containsKey(arr[1])){            		            		
            	
            		genes.put(arr[1],"BioCyc");
            		
            		phge.put(arr[3],genes);
            	
            	}                        	
            	            	
            }else{
            	
            	genes.put(arr[1],"BioCyc");
            	
            	phge.put(arr[3],genes);
            }
            
            
            Pattern idp = Pattern.compile("MGI:(\\d{3,15})");
            
            Matcher idm = idp.matcher(ln);
            
            if(idm.find()){
            	
            	Vector<String> pv = new Vector<String>();
            	
                id = idm.group(1);
                //System.out.println(id);            
            
                if(!idmp.containsKey(id)){
                	
                    idmp.put(id, arr[1]);
                    
                }
                
                arr[3] = arr[3].trim();
                
                if(!pathmp.containsKey(id)){  
                	
                    pv.add(arr[3]);
                    
                    pathmp.put(id,pv);
                    
                }else{
                    
                	pv = pathmp.get(id);
                	
                    if(!pv.contains(arr[3])){
                    	
                        pv.add(arr[3]);
                        
                    }
                    
                    pathmp.put(id,pv);
                }
            }
        }
        br.close();
        System.out.println( "The size of node ID map is: "+idmp.size()+"\n"+ "The size of node pathway name map is: "+pathmp.size());
        
    }

    private static void Writepairrelwithpathname(HashMap<String, String> idmp, HashMap<String, Vector> pathmp, String pair, String res)throws IOException {
        
        String ln;
        String[] arr;        
        
        BufferedReader br = new BufferedReader(new FileReader(pair));
        PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(res)));
        
        HashMap<String,Boolean> ckmp = new HashMap<String,Boolean>();

        while((ln = br.readLine())!=null){
            
        	arr = ln.split("\t");
            
            
            if(pathmp.containsKey(arr[6]) && pathmp.containsKey(arr[7])){
                
                //System.out.println("debug");
                Vector pva = pathmp.get(arr[6]);
                
                Vector pvb = pathmp.get(arr[7]);
                
                Vector<String> ckv = new Vector<String>();
                
                for(int a = 0;a<pva.size();a++){
                    
                	for(int b = 0;b <pvb.size();b++){
                    
                		if(pva.get(a).equals(pvb.get(b))){
                        
                			if(!ckv.contains(((String) pva.get(a)).trim())){
                                
                            	ckv.add(((String) pva.get(a)).trim());
                                
                                String geneA = idmp.get(arr[6]);
                                
                                String geneB = idmp.get(arr[7]);
                                //if(!ckmp.containsKey(geneA+"\t"+geneB)&&!ckmp.containsKey(geneB+"\t"+geneA)){
                                	
                                	String pathway = (String) pva.get(a);
                                    
                                	pw.println(geneA+"\t"+geneB+"\t"+arr[5]+"\t"+pathway);
                                    
                                	ckmp.put(geneA+"\t"+geneB, Boolean.TRUE);                                                                        
                                    
                                    HashMap<String,String> pairs = new HashMap<String,String>();
                                    
                                    String ga, gb;
                                    
                                    if(geneA.compareTo(geneB)>=0){
                                    	
                                    	ga = geneA;
                                    	
                                    	gb = geneB;
                                    	
                                    }else{
                                    	
                                    	gb = geneA;
                                    	
                                    	ga = geneB;
                                    	
                                    }
                                    
                                    if(phgp.containsKey(pva.get(a))){
                                    	
                                    	pairs = phgp.get(pathway);
                                    	
                                    	if(!pairs.containsKey(ga+"\t"+gb)){
                                    		
                                    		pairs.put(ga+"\t"+gb, arr[5]);
                                    		
                                        	phgp.put(pathway, pairs);
                                    	}                                    	
                                    	
                                    }else{
                                    	
                                    	pairs.put(ga+"\t"+gb, arr[5]);
                                    	
                                    	phgp.put(pathway, pairs);
                                    	
                                    }
                                //}
                            }                           
                        }
                    }
                }                
            }
        }
        
        br.close();
        
        pw.close();
    }

    
}