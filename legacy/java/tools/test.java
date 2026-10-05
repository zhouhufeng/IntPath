package tools;

import java.util.*;
import java.io.*;

public class test {
	
	public void writeEveryThingOut(HashMap<String, HashMap<String, String>> mp)throws Exception{
		
		Set kset = mp.keySet();
		
		for(Object objk : kset){
			
			String k = objk.toString();
			
			HashMap<String,String> vmp = new HashMap<String,String>();
			
			vmp = mp.get(k);
			
			Set vset = vmp.keySet();
			
			for(Object objvk: vset){
				
				String vk = objvk.toString();
				
				String vv = vmp.get(vk);
				
				System.out.println(k+"\t"+vk+"\t"+vv);
				
			}
			
		}
		
	}
	
	public void writeOutVectorKey( HashMap<String, Vector<String>>mp)throws Exception{
		
		Set kset = mp.keySet();
		
		for(Object objk : kset){
			
			String k = objk.toString();
			
			Vector<String> vmp = new Vector<String>();
			
			vmp = mp.get(k);
			
			for(Object objvk: vmp){
				
				String vk = objvk.toString();
				
				System.out.println(k+"\t"+vk);
				
			}
			
		}
		
	}	
	
	public void writempvToFile( HashMap<String, Vector<String>>mp,String f)throws Exception{
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f)));
		
		Set<String> kset = mp.keySet();
		
		for(String k : kset){			
			
			Vector<String> vmp = new Vector<String>();
			
			vmp = mp.get(k);
			
			//for(String vk : vmp){		
		    for(int i = 0; i < vmp.size(); i++){
				
				String vk = vmp.get(i);
				//System.out.println(k+"\t"+vk);
				
				pw.println(k+"\t"+vk);
			}
			
		}
		
		pw.close();
	}
	
	public void writeGenesToFile(HashMap<String, HashMap<String, String>> mp, String otfl)throws IOException{
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(otfl)));
		
		Set kset = mp.keySet();
		
		for(Object objk : kset){
			
			String pwy = objk.toString();
			
			HashMap<String,String> vmp = new HashMap<String,String>();
			
			vmp = mp.get(pwy);
			
			Set vset = vmp.keySet();
			
			for(Object objvk : vset){
				
				String gen = objvk.toString();
				
				String db = vmp.get(gen);
				
				//System.out.println(gen+"\t"+pwy+"\t"+db);
				
				//pw.println(gen+"\t"+pwy+"\t"+db);
				
				pw.println(pwy+"\t"+gen+"\t"+db);
				
			}
			
		}
		
		pw.close();
		
	}	
		
	public void writeGenePairsToFile(HashMap<String, HashMap<String, String>> mp, String otfl, String db)throws IOException{
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(otfl)));
		
		Set kset = mp.keySet();
		
		for(Object objk : kset){
			
			String pwy = objk.toString();
			
			HashMap<String,String> vmp = new HashMap<String,String>();
			
			vmp = mp.get(pwy);
			
			Set vset = vmp.keySet();
			
			for(Object objvk : vset){
				
				String genpair = objvk.toString();
				
				String rel = vmp.get(genpair);
				
				//System.out.println(genpair+"\t"+rel+"\t"+pwy+"\t"+db);
				
				pw.println(genpair+"\t"+rel+"\t"+pwy+"\t"+db);
				
			}
			
		}
		
		pw.close();
		
	}	
	
	public static void writeRelPWYgroup(HashMap<String,Vector<String>> mp, PrintWriter pw)throws IOException{
		
		//PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
		
		Set mpkset = mp.keySet();
		
		int pid = 0;
		
		for(Object kO:mpkset){
			
			String k = kO.toString();
			
			Vector<String> v = new Vector<String>();
			
			v= (Vector<String>) mp.get(k);					
			
			for(int i = 0; i< v.size(); i++){
				
				pid ++;
				
				String vs = v.get(i);
				
				pw.println(pid+"\t"+k+"\t"+vs);
				
				//System.out.println(pid+"\t"+k+"\t"+vs);
				
			}
			
		}
		
		//pw.close();
		
		
	}
	
	public int NumHashInHash(HashMap<String, HashMap<String, String>> mp) {
		
		int tnum = 0;
		
		Set<String> kset = mp.keySet();		
		
		for(String k : kset){
			
			HashMap<String,String> vmp = mp.get(k);
			
			tnum = tnum + vmp.size();
			
		}
		
		return tnum;
		
	}

	public void writeGenePairsToFile( HashMap<String, HashMap<String, String>> mp, String otfl) throws IOException {
		
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(otfl)));
		
		Set<String> kset = mp.keySet();
		
		for(String pwy : kset){			
			
			HashMap<String,String> vmp = new HashMap<String,String>();
			
			vmp = mp.get(pwy);
			
			Set<String> vset = vmp.keySet();
			
			for(String gep : vset){
				
				String reldb = vmp.get(gep);
				
				String[] rd = reldb.split("\t");
				
				pw.println(gep+"\t"+rd[0]+"\t"+pwy+"\t"+rd[1]);
				
				//pw.println(pwy+"\t"+gep+"\t"+db);
				
			}
			
		}
		
		pw.close();
		
	}

	public static void OutputGroups(HashMap<String, HashMap<String, Vector<String>>> gro,
			String db, PrintWriter pw) {
		/*
		 *  Output the genes information in the group,
		 *  format, Pathways, genes, groupID, source database
		 */
		Set<String> pwyset = gro.keySet();
		
		for(String pwy: pwyset){
			
			HashMap<String,Vector<String>> vmp = gro.get(pwy);
			
			Set<String> groidset = vmp.keySet();
			
			for(String groid:groidset){
				
				Vector<String> vv = vmp.get(groid);
				
				for(int i = 0; i<vv.size();i++){
					
					String ge = vv.get(i);
					
					pw.println(pwy+"\t"+ge+"\t"+groid+"\t"+db);
					
				}
				
			}
			
		}
	
		
	}

	public static void MappingContainers(HashMap<String, Vector<String>> ecmp) {
		// This method is used to test the mapping files.
		
		Set<String> otherset = ecmp.keySet();
		
		for(String others: otherset){
			
			Vector<String> tvv = ecmp.get(others);
			
			if(tvv.size()>1){
				
				for(String ts: tvv){
					
					System.out.print(ts);
					
				}
				
				System.out.print("\n");
				
			}
			
		}
		
	}
	
}
