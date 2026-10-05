package utils;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.*;

public class stats {
	
	public stats(){
		
	}
	
	public int calcGenes(String organism) {
		
		Vector<String> qus = new Vector<String>();
		
		qus.add("gene_names");
		
		String sqlQuery = "select gene_names from pathway_genes";
		
		return query(organism, qus, sqlQuery);
	}
	
	public int calcGenePairs(String organism) {

		Vector<String> qus = new Vector<String>();
		
		String sqlQuery = "select gene_A,gene_B from gene_pairs";
		
		qus.add("gene_A");
		
		qus.add("gene_B");
				
		return query(organism, qus, sqlQuery);
	}
	
	public int calcIntPath(String organism){
		
		Vector<String> qus = new Vector<String>();
		
		qus.add("integrated_pathway_name");
		
		String sqlQuery = "select integrated_pathway_name from related_pathways";		

		return query(organism, qus, sqlQuery);
	}
	
	public int calcOrigiPath(String organism){
		
		Vector<String> qus = new Vector<String>();
		
		qus.add("pathway_id");
		
		String sqlQuery = "select pathway_id from pathway_names";		

		return query(organism, qus, sqlQuery);
	}	

	
	public int calcKEGG (String organism){
		
		Vector<String> qus = new Vector<String>();
		
		qus.add("pathway_name");
		
		String sqlQuery = "select pathway_name,pathway_origin from pathway_names where pathway_names.pathway_origin = 'KEGG' ";		

		return query(organism, qus, sqlQuery);
	}		
	
	public int calcWiki (String organism){
		
		Vector<String> qus = new Vector<String>();
		
		qus.add("pathway_name");
		
		String sqlQuery = "select pathway_name,pathway_origin from pathway_names where pathway_names.pathway_origin = 'WikiPathways' ";		

		return query(organism, qus, sqlQuery);
	}	
	
	public int calcBioCyc (String organism){
		
		Vector<String> qus = new Vector<String>();
		
		qus.add("pathway_name");
		
		String sqlQuery = "select pathway_name,pathway_origin from pathway_names where pathway_names.pathway_origin = 'BioCyc' ";		

		return query(organism, qus, sqlQuery);
	}		
	
	public int query (String organism, Vector<String> qus, String sqlQuery){
	
		int c = 0;
		
		Connection con = null;
		Statement smtm = null;
		ResultSet rs   = null;
	
		try {
			
			Class.forName("com.mysql.jdbc.Driver"); 
			
			if (organism.equalsIgnoreCase("M.musculus")){
				con = DriverManager.getConnection("jdbc:mysql://localhost:3666/musculus?user=INTPATH_DB_USER&password=INTPATH_DB_PASSWORD");
			}
			else if (organism.equalsIgnoreCase("H.sapiens")){
				con = DriverManager.getConnection("jdbc:mysql://localhost:3666/sapiens?user=INTPATH_DB_USER&password=INTPATH_DB_PASSWORD");
			}
			else if (organism.equalsIgnoreCase("S.cerevisiae")){
				con = DriverManager.getConnection("jdbc:mysql://localhost:3666/cerevisiae?user=INTPATH_DB_USER&password=INTPATH_DB_PASSWORD");
			}
			else if (organism.equalsIgnoreCase("M.tuberculosis_H37Rv")){
				con = DriverManager.getConnection("jdbc:mysql://localhost:3666/tuberculosis?user=INTPATH_DB_USER&password=INTPATH_DB_PASSWORD");
			}		
	
			smtm = con.createStatement();
			
			rs   = smtm.executeQuery(sqlQuery);
			
			HashMap<String,Integer> ctmp = new HashMap<String,Integer>();
			
			int idx = 0;
			
			while (rs != null && rs.next()) {
				
				String ot = "";
				
				for(int i = 0;i< qus.size();i++){
					
					String tmp = rs.getString(qus.get(i));
					
					ot = ot +"\t"+ tmp;
					
				}				
				
				idx++;
				
				ctmp.put(ot, idx);									
				
			}
			
			c = ctmp.size();
			
		}catch (ClassNotFoundException e){
			
			e.printStackTrace();
			
		}catch (SQLException ex) {
			
			System.out.println(ex.getMessage());
			System.out.println(ex.getSQLState());
			System.out.println(ex.getErrorCode());
	
		} finally {
			
			try {
				
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
				
			} catch (Exception ee) {
				
			}
		}
	
		return c;
		
	   }	
}
