package tools;

import java.util.*;


public class calc {
	
	public void calc(){				
	}
	/*
	 * Note: this method is specifically used to handle the different number of pathways in pathway genes 
	 * and pathway-gene pairs containers. 
	 */
	
	protected double avge, avgp, ngp, nge, npge, npgp;
	
	public void  AverGENnGPR(HashMap<String, HashMap<String, String>> ge, 
			HashMap<String, HashMap<String, String>> gp) {
		
		double avge, avgp, ngp, nge, npge, npgp;
		
		int tge, tgp;
		
		tge = tgp = 0;
		
		avgp = avge = avgp = ngp = nge = npge = npgp = 0.0;
		
		npge = ge.size(); //number of pathways in pathway genes
		
		npgp = gp.size(); //number of pathways in pathway gene-pairs
				
		
		//----iterate gene pairs -------
		Set<String> phgpset = gp.keySet();	
		
		for(String pwygp : phgpset){		
			
			int pwygpn = gp.get(pwygp).size();
			
			tgp = tgp + pwygpn;			
			
		}		
		
		avgp = tgp / npge;
		
		//-----iterate genes-----------
		
		Set<String> phgeset = ge.keySet();	
		
		for(String pwyge : phgeset){		
			
			int pwygen = ge.get(pwyge).size();
			
			tge = tge + pwygen;			
			
		}		
		
		avge = tge / npge;
		
		System.out.println(
				"Specifically, total pairs number: "+tgp+
				". \n total genes number: "+tge+
				". \n All the pathways number is(pathway-genes): "+npge+
				". \n The number of pathways in pathway-gene pairs: "+npgp+
				". \n Average number of gene per pathway: "+avge+
				". \n Average number of gene pairs per pathway: "+avgp+".\n"
				);
	}
	
	public double AverNUMperPWY(HashMap<String, HashMap<String, String>> pg) {
				
		Set phset = pg.keySet();	
		
		double aver = 0.00;
		
		int total = 0;
		
		for(Object pathwaygens: phset){
			
			String pwy = pathwaygens.toString();
			
			int pwyg = pg.get(pwy).size();
			
			total = total + pwyg;			
			
		}		
		
		aver = total / pg.size();
		
		System.out.println(
				"Specifically, total genes/pairs: "+total+
				". All the pathways number is: "+pg.size()
				);
		// TODO Auto-generated method stub
		
		return aver;
	}			
    
    public static void compGenes(HashMap amp,String ta, HashMap bmp, String tb){
        
        Set aks = amp.keySet();
        
        int a = amp.size();
        
        int b = bmp.size();
        
        int o = 0;
        
        for(Object ao:aks){
        
        	String ak = ao.toString();
            
        	if(bmp.containsKey(ak)){
            
        		o++;
            
        	}
        
        }
        
        System.out.println(
                "The number of genes in "+ta+" : "+a+
                "\nThe number of genes in "+tb+" : "+b+
                "\nThe number of overlap genes between "+ta+" and "+tb+"is: "+o+
                "\nThe number of unique genes between "+ta+" and "+tb+"is: "+(a+b-2*o)+
                "\nJaccard coefficient between "+ta+" and "+tb+"is: "+ ((double)o/(a+b-o))
                );        
    }
    
    public static void compPairs(HashMap amp,String ta, HashMap bmp, String tb){
        
        int a = amp.size();
        
        int b = bmp.size();
        
        int o = 0;        
        
        Set aset = amp.keySet();
        
        for(Object aso:aset){
            
        	String pair = aso.toString();
            
            String[] p  = pair.split("\t");
            
            if(bmp.containsKey(p[0]+"\t"+p[1])||bmp.containsKey(p[1]+"\t"+p[0])){
            
            	o++;
            	
            }
        }
        
        System.out.println(
                "The number of gene pairs in "+ta+" :"+a+
                "\nThe number of gene pairs in "+tb+" :"+b+
                "\nThe number of overlap gene pairs between "+ta+" and "+tb+"is: "+o+
                "\nThe number of unique gene pairs between "+ta+" and "+tb+"is: "+(a+b-2*o)+
                "\nJaccard coefficient between "+ta+" and "+tb+"is: "+ ((double)o/(a+b-o))
                );             
    }


    public static void compDBGenes(	HashMap<String, HashMap<String, String>> kge, 
    		String dbk,	HashMap<String, HashMap<String, String>> cge, String dbc) {
    	
		HashMap<String,String> kmp = new HashMap<String,String>();
		
		HashMap<String,String> cmp = new HashMap<String,String>();
		
		Set<String> kset = kge.keySet();
		
		for(String ks : kset){
			
			kmp.putAll(kge.get(ks));
			
		}
		
		
		Set<String> cset = cge.keySet();
		
		for(String cs : cset){
			
			cmp.putAll(cge.get(cs));
			
		}
		
		compGenes(kmp,dbk,cmp,dbc);
	}

	public static void compDBGenePairs(HashMap<String, HashMap<String, String>> kgp, 
			String dbk,HashMap<String, HashMap<String, String>> cgp, String dbc) {
		// TODO Auto-generated method stub
		
		HashMap<String,String> kmp = new HashMap<String,String>();
		
		HashMap<String,String> cmp = new HashMap<String,String>();
		
		Set<String> kset = kgp.keySet();
		
		for(String ks : kset){
			
			kmp.putAll(kgp.get(ks));
			
		}
		
		Set<String> cset = cgp.keySet();
		
		for(String cs : cset){
			
			cmp.putAll(cgp.get(cs));
			
		}
		
		compPairs(kmp,dbk,cmp,dbc);
		
	}

	public int threeOverlapPWYs(HashMap<String, Vector<String>> intpwyname) {
		
		/*
		 *  By calling this method will return the number of pathway group have 
		 *  pathways from three source database.
		 */
		
		int n = 0;
		
		Set<String> intpathnameset = intpwyname.keySet();
		
		for(String keynames: intpathnameset){
			
			Vector<String> orinmsv = intpwyname.get(keynames);
			
			int k , c, w;
			k = c = w = 0;		
			
			for(String os : orinmsv){
				
				String[] tag = os.split("\\+");
				
				if(tag[1].equalsIgnoreCase("k")){
					
					k++;
					
				}else if(tag[1].equalsIgnoreCase("C")){
					
					c++;
					
				}else if(tag[1].equalsIgnoreCase("W")){
					
					w++;
					
				}
				
				
			}
			
			if(k>0 && c>0 && w>0){
				
				n++;
				
			}
			
		}
		
		return n;
	}


}
