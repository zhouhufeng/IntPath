package tools;

import java.util.*;
import java.io.*;

public class conts {
	
	/*
	 * This containers is used for the Integration class.
	 * input: HashMap<String，HashMap<String,String>> pathwayGenes and 
	 * 		  HashMap<String，HashMap<String,String>> pathwayGenePaires
	 * 
	 * output: Vector<String> pathwayNamesGenes and
	 * 		   Vector<String> pathwayNamesGenePairs
	 * 
	 */
	
	public static void main(String args[])throws IOException{
		
		preparePPIN();
		
	}
	
	private static void preparePPIN() throws IOException {
		// This method is used to prepare the PPIN used in calculating the distances between pathways.
		
		String orgs = "musculus";
		
		String outmmu = orgs+File.separator+"source"+File.separator+"PPI"+File.separator+orgs+"STRING";		
		
		int thrldmmu = 750;
		
		isoSTRINGPPI(orgs, "10090.", outmmu,thrldmmu);
			
		
		
		orgs = "sapiens";
				
		String outhsa = orgs+File.separator+"source"+File.separator+"PPI"+File.separator+orgs+"STRING";
		
		int thrldhsa = 750;
		
		isoSTRINGPPI(orgs, "9606.", outhsa,thrldhsa);
						
		
		orgs = "cerevisiae";
		
		
		String outces = orgs+File.separator+"source"+File.separator+"PPI"+File.separator+orgs+"STRING";
		
		int thrldces = 750;
		
		isoSTRINGPPI(orgs, "4932.", outces,thrldces);
		
		
		orgs = "tuberculosis";	
		
		
		String outmtb = orgs+File.separator+"source"+File.separator+"PPI"+File.separator+orgs+"STRING";
		
		int thrldmtb = 750;
		
		isoSTRINGPPI(orgs, "83332.", outmtb,thrldmtb);
		
	}
	
	
	public void pathwayContainers(){
		
		
	}
	
	
	
	public Vector<String> keyToVec(HashMap<String,HashMap<String,String>> mp){
		
		Vector<String> vec = new Vector<String>();
		
		Set mpkeyset = mp.keySet();
		
		for(Object keyO : mpkeyset){
			
			String mpkeys = keyO.toString();
			
			vec.add(mpkeys);
			
		}
		
		return vec;
	}
	
	public static void isoSTRINGPPI (String orgs, String orgID,String outp,int thrld)throws IOException{
		
		String inputf     = "rawdata"+File.separator+"PPI"+File.separator+"STRINGPPIsuptOrgs";
		
		BufferedReader br = new BufferedReader(new FileReader(inputf));
		
		PrintWriter pw    = new PrintWriter(new BufferedWriter (new FileWriter(outp)));
		
		String ln;
		
		String[] arr;
		
		HashMap<String,String> ck = new HashMap<String,String>();
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split(" ");
			
			int fs = Integer.parseInt(arr[15]);
			
			String[] ga = arr[0].split("\\.");
			
			if(ga.length!=2)continue;
									
			String[] gb = arr[1].split("\\.");
			
			if(gb.length!=2)continue;
			
			if(ln.contains(orgID) && fs>=700){
				
				//pw.println(ln);
				
				String gena,genb;
				
				if(ga[1].compareTo(gb[1])>=0){
					
					gena = ga[1];
					
					genb = gb[1];
					
				}else{
					
					genb = ga[1];
					
					gena = gb[1];
					
				}
				
				if(!ck.containsKey(gena+"\t"+genb)){
					
					ck.put(gena+"\t"+genb, null);
					
					pw.println(gena+"\t"+genb);
					
				}
				
			}
		}
		
		
		br.close();
		
		pw.close();
	}

	public static HashMap<String,String> PPImaping(HashMap<String, Vector<String>> idmp,
			String outmmu, String cleanPPI) throws IOException{
		
		String ln;
		
		String arr[];//,ga,gb;
		
		HashMap<String,String> ppi = new HashMap<String,String>();
		
		BufferedReader br = new BufferedReader(new FileReader(outmmu));
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(cleanPPI)));
		
		while((ln = br.readLine())!=null){
			
			arr = ln.split("\t");
			
//			String geneA = arr[0];
//			
//			String geneB = arr[1];
//			
//			String[] genea = geneA.split("\\.");
//			
//			String[] geneb = geneB.split("\\.");
			
			if(idmp.containsKey(arr[0]) && idmp.containsKey(arr[1])){
				
				Vector<String> va = idmp.get(arr[0]);
				
				Vector<String> vb = idmp.get(arr[1]);
				
				for(int i = 0; i<va.size();i++){
					
					for(int j = 0; j<vb.size();j++){
						
						String ga = va.get(i);
						
						String gb = vb.get(j);
						
						String a = "";
						
						String b = "";
						
						if(ga.compareTo(gb)>=0){
							
							a = ga;
							
							b = gb;
							
						}else{
							
							a = gb;
							
							b = ga;
														
						}
						
						if (a == ""||b == "") continue;
						
						if(!ppi.containsKey(a+"\t"+b)){
							
							ppi.put(a+"\t"+b, null);
							
							pw.println(a+"\t"+b);
							
						}
						
					}
					
				}
				
			}
			
			
		}
		
		br.close();
		
		pw.close();
		
		return ppi;
		
	}

}
