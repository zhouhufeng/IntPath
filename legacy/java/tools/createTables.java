package tools;

import java.io.*;
import java.sql.*;
import java.util.*;

public class createTables {
	
	static Connection con  = null; 
	static Statement  smtm = null;
	static ResultSet  rs   = null;	
	
	
//	public void createTables(){
//		
//		
//		Class.forName("com.mysql.jdbc.Driver"); 
//		
//		con  = DriverManager.getConnection("jdbc:mysql://localhost/"+org+"?user=INTPATH_DB_USER&password=INTPATH_DB_PASSWORD");
//		    
//		smtm = con.createStatement();	
//		
//		
//		
//		if (rs != null) {
//			rs.close();
//			rs = null;
//		}
//		if (smtm != null) {
//			smtm.close();
//			smtm = null;
//		}
//		if (con != null) {
//			con.close();
//			con = null;
//		}		
//				
//		
//	}
	
	/*

	 */
	
	public static void ctGeneMapping(String org,HashMap<String,Vector<String>> mp) 
			throws IOException, ClassNotFoundException, SQLException{
		
		Class.forName("com.mysql.jdbc.Driver"); 
		
		con  = DriverManager.getConnection("jdbc:mysql://localhost/"+org+"?user=INTPATH_DB_USER&password=INTPATH_DB_PASSWORD");
		    
		smtm = con.createStatement();	
		
		smtm.executeUpdate("delete from gene_mapping");	
		
		int gid = 0;
		
		if(mp == null) return;
		
		Set onameset = mp.keySet(); //oname stands for other names, kname stands for standard name.
		
		for(Object onameO : onameset){
			
			String oname = onameO.toString();						
			
			Vector<String> kv = null;						
			
			kv = mp.get(oname);		
			
			if(kv == null) continue;
			
			if(kv.isEmpty())continue;	
			
			oname = oname.replace("'", "\\'"); //lnh = lnh.replace("'", "\\'");
			
			for(int i = 0; i<kv.size(); i++){
				
				String kname = kv.get(i);
				
				if(kname.equals("?")||oname.equals("?")){
					continue;
				}
				
				//System.out.println(oname);
				
				if(oname.length()>80){
					
					System.out.println("This name won't be insert into gene mapping Mysql table: "+ oname);
					
					continue;
					
				}
				
				gid++;
				
				String ctgm = "INSERT INTO gene_mapping (gene_id , unified_gene_name , other_gene_names ) VALUES (" + gid + ", '" + kname + "', '"+ oname+ "')";
				
				//System.out.println(gid + "\t" + kname+ "\t" +oname);

				smtm.executeUpdate(ctgm);	
				
			}
			
			
		}
		
		if (rs != null) {
			rs.close();
			rs = null;
		}
		if (smtm != null) {
			smtm.close();
			smtm = null;
		}
		if (con != null) {
			con.close();
			con = null;
		}			

	}


	public static void ctPathwayNames(Vector<String> cv,Vector<String> kv, Vector<String> wv, 
			String orgs) throws IOException, ClassNotFoundException, SQLException {
		
		Class.forName("com.mysql.jdbc.Driver"); 
		
		con  = DriverManager.getConnection("jdbc:mysql://localhost/"+orgs+"?user=INTPATH_DB_USER&password=INTPATH_DB_PASSWORD");
		    
		smtm = con.createStatement();	
		
        smtm.executeUpdate("delete from pathway_names");
        
        int id = 0;        
		
		for(int c = 0; c<cv.size(); c++){
			
			String pwy = cv.get(c);
			
			String sdb = "BioCyc";
			
			if(pwy.contains("'")){
				
				pwy = pwy.replace("'", "\\'"); 
				
			}else if(pwy.contains("\'")){
				
				pwy = pwy.replace("\'", "\\'"); 
				
			}
			
			id++;
			
			smtm.executeUpdate("INSERT INTO pathway_names(pathway_id, pathway_name, pathway_origin) VALUES("+id+",'"+pwy+"','"+sdb+"')");
		}
		
		for(int k = 0; k<kv.size(); k++){
			
			String pwy = kv.get(k);
			
			String sdb = "KEGG";
			
			if(pwy.contains("'")){
				
				pwy = pwy.replace("'", "\\'"); 
				
			}else if(pwy.contains("\'")){
				
				pwy = pwy.replace("\'", "\\'"); 
				
			}			
			
			id++;
			
			smtm.executeUpdate("INSERT INTO pathway_names(pathway_id, pathway_name, pathway_origin) VALUES("+id+",'"+pwy+"','"+sdb+"')");
		}
		
		for(int w = 0; w<wv.size(); w++){
			
			String pwy = wv.get(w);
			
			String sdb = "WikiPathways";
			
			if(pwy.contains("'")){
				
				pwy = pwy.replace("'", "\\'"); 
				
			}else if(pwy.contains("\'")){
				
				pwy = pwy.replace("\'", "\\'"); 
				
			}
			
			id++;
			
			smtm.executeUpdate("INSERT INTO pathway_names(pathway_id, pathway_name, pathway_origin) VALUES("+id+",'"+pwy+"','"+sdb+"')");
			
		}
		
		
		if (rs != null) {
			rs.close();
			rs = null;
		}
		if (smtm != null) {
			smtm.close();
			smtm = null;
		}
		if (con != null) {
			con.close();
			con = null;
		}	
		
	}


	public static void ctPathwayGenes(HashMap<String, HashMap<String, String>> intPathGEN, 
			String orgs) throws IOException, ClassNotFoundException, SQLException {
		
		Class.forName("com.mysql.jdbc.Driver"); 
		
		con  = DriverManager.getConnection("jdbc:mysql://localhost/"+orgs+"?user=INTPATH_DB_USER&password=INTPATH_DB_PASSWORD");
		    
		smtm = con.createStatement();	
		
		smtm.executeUpdate("delete from pathway_genes");
		
		
		Set<String> kset = intPathGEN.keySet();
		
		for(String pwy : kset){
			
			
			HashMap<String,String> vmp = new HashMap<String,String>();
			
			vmp = intPathGEN.get(pwy);
			
			Set<String> vset = vmp.keySet();
			
			for(String gen: vset){
				
				
				String db = vmp.get(gen);
				
				if(pwy.contains("\\'")){
					
				}else if(pwy.contains("'")){
					
					pwy = pwy.replace("'", "\\'"); 
					
				}else if(pwy.contains("\'")){
					
					pwy = pwy.replace("\'", "\\'"); 
					
				}								
				
				
				//System.out.println(pwy+"\t"+gen+"\t"+db);
				
				String ctpthgen = "INSERT INTO pathway_genes (pathway_name , gene_names , origin_database ) VALUES ('" + pwy + "', '" + gen + "', '"+ db+ "')";
				
				smtm.executeUpdate(ctpthgen);								
				
			}
			
		}		
		
		
		
		if (rs != null) {
			rs.close();
			rs = null;
		}
		if (smtm != null) {
			smtm.close();
			smtm = null;
		}
		if (con != null) {
			con.close();
			con = null;
		}	
		
	}


	public static void ctPathwayGenePairs(HashMap<String, HashMap<String, String>> mp, 
			String orgs) throws IOException, ClassNotFoundException, SQLException {

		Class.forName("com.mysql.jdbc.Driver"); 
		
		con  = DriverManager.getConnection("jdbc:mysql://localhost/"+orgs+"?user=INTPATH_DB_USER&password=INTPATH_DB_PASSWORD");
		    
		smtm = con.createStatement();	
		
		smtm.executeUpdate("delete from pathway_allpairinfo");
		
		smtm.executeUpdate("delete from gene_pairs");
				
		
		Set<String> kset = mp.keySet();
		
		for(String pwy : kset){			
			
			HashMap<String,String> vmp = new HashMap<String,String>();
			
			vmp = mp.get(pwy);
			
			Set<String> vset = vmp.keySet();
			
			for(String gep : vset){
				
				String reldb = vmp.get(gep);
				
				String[] rd = reldb.split("\t");
				
				String[] ge = gep.split("\t");
				
				//System.out.println(gep+"\t"+rd[0]+"\t"+pwy+"\t"+rd[1]);
				
				
				if(pwy.contains("\\'")){
					
				}else if(pwy.contains("'")){
					
					pwy = pwy.replace("'", "\\'"); 
					
				}else if(pwy.contains("\'")){
					
					pwy = pwy.replace("\'", "\\'"); 
					
				}								
				
				
				String s1 = "INSERT INTO pathway_allpairinfo (gene_A , gene_B , relations , pathway_name , source_database) VALUES ('" + ge[0] + "', '" + ge[1] + "', '" + rd[0] + "', '" + pwy + "', '" + rd[1] + "')";
				//System.out.println(s);
				
				smtm.executeUpdate(s1);
				
				String s2 = "INSERT INTO gene_pairs (gene_A , gene_B , relations , pathway_name ) VALUES ('" + ge[0] + "', '" + ge[1] + "', '" + rd[0] + "', '" + pwy + "')";
				
				smtm.executeUpdate(s2);							
				
			}
			
		}		
		
		
		
		if (rs != null) {
			rs.close();
			rs = null;
		}
		if (smtm != null) {
			smtm.close();
			smtm = null;
		}
		if (con != null) {
			con.close();
			con = null;
		}
		
	}


	public static void ctRelatedPathways( HashMap<String, Vector<String>> relPth, 
			String orgs)throws IOException, ClassNotFoundException, SQLException {
		
		Class.forName("com.mysql.jdbc.Driver"); 
		
		con  = DriverManager.getConnection("jdbc:mysql://localhost/"+orgs+"?user=INTPATH_DB_USER&password=INTPATH_DB_PASSWORD");
		    
		smtm = con.createStatement();	
		
		smtm.executeUpdate("delete from related_pathways");
		
		Set<String> intpwyset = relPth.keySet();
		
		for(String intpwy : intpwyset ){
			
			Vector<String> pwyv = relPth.get(intpwy);			
			
			HashMap<String,String> intgens = new HashMap<String,String>();	
			
			int ctrmv = 0;
			
			for(int i = 0; i < pwyv.size(); i++){
				
				String pwydb = pwyv.get(i);
				
				String[] pdb = pwydb.split("\\+");
				
				//HashMap<String,HashMap<String,String>> ge = null;
							
				HashMap<String,String> epge = null;
				
				String pwy = pdb[0];
				
				if(pwy.contains("'")){
					
					pwy = pwy.replace("'", "\\'"); 
					
				}else if(pwy.contains("\'")){
					
					pwy = pwy.replace("\'", "\\'"); 
					
				}
				
				if(pdb[1].equalsIgnoreCase("K")){															
					
					smtm.executeUpdate("INSERT INTO related_pathways(integrated_pathway_id, integrated_pathway_name, pathway_name, source_database) VALUES("+ctrmv+",'"+intpwy+"','"+pwy+"','KEGG')");
					
				}else if (pdb[1].equalsIgnoreCase("W")){
				
					smtm.executeUpdate("INSERT INTO related_pathways(integrated_pathway_id, integrated_pathway_name, pathway_name, source_database) VALUES("+ctrmv+",'"+intpwy+"','"+pwy+"','WikiPathways')");                    
		            
				}else if (pdb[1].equalsIgnoreCase("C")){
					
					smtm.executeUpdate("INSERT INTO related_pathways(integrated_pathway_id, integrated_pathway_name, pathway_name, source_database) VALUES("+ctrmv+",'"+intpwy+"','"+pwy+"','BioCyc')");					  
					
				}else{
					
					System.out.println("Warning! In BuildIntPathGEN Some original pthway names do not have + tags");
					
				}				
			}
		}
		
		if (rs != null) {
			rs.close();
			rs = null;
		}
		if (smtm != null) {
			smtm.close();
			smtm = null;
		}
		if (con != null) {
			con.close();
			con = null;
		}

	}
}