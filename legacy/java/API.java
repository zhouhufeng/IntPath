package API;

import java.io.*;
import java.sql.*;
import java.util.*;

public class API {

	public String getGeneID(String name,String organism) {
		String typeName = "unified_gene_name";
		String sqlQuery = "select unified_gene_name from gene_mapping where gene_mapping.other_gene_names=\""+ name + "\"";
		return getResult(organism, typeName, sqlQuery);
	}

	public String getDBPathways(String databasename,String organism) {
		//System.out.println("debug");
		List<String> typeIDList = new ArrayList<String>();
		
		typeIDList.add("DatabaseName: ");
		typeIDList.add("PathwayID: ");
		typeIDList.add("PathwayName: ");
		
		List<String> typeNameList = new ArrayList<String>();
		
		typeNameList.add("pathway_origin");
		typeNameList.add("pathway_id");
		typeNameList.add("pathway_name");
		
		String sqlQuery = "select pathway_id,pathway_name,pathway_origin from pathway_names where pathway_names.pathway_origin =\""+databasename+"\"";
		
		return getResult(organism, typeIDList, typeNameList, sqlQuery);
	}

	public String getPathway(String pathwayname,String organism) {
		
		List<String> typeIDList = new ArrayList<String>();
		typeIDList.add("DatabaseName: ");
		typeIDList.add("PathwayID: ");
		typeIDList.add("PathwayName: ");
		
		List<String> typeNameList = new ArrayList<String>();
		typeNameList.add("pathway_origin");
		typeNameList.add("pathway_id");
		typeNameList.add("pathway_name");
		
		String sqlQuery = "select pathway_id,pathway_name,pathway_origin from pathway_names where pathway_names.pathway_name=\""+pathwayname+"\"";
			
		return getResult(organism, typeIDList, typeNameList, sqlQuery);
		
	}

	public String getPathwayGene(String pathwayname,String organism) {
		String typeName = "gene_names"; 
		String sqlQuery = "select gene_names from pathway_genes where pathway_genes.pathway_name=\""+ pathwayname + "\"";
		return getResult(organism, typeName, sqlQuery);
	}

	public String getGenePathway(String geneID,String organism) {
		String sqlQuery = "select pathway_name from pathway_genes where pathway_genes.gene_names=\""+ geneID + "\"";
		String typeName = "pathway_name";
		return getResult(organism, typeName, sqlQuery);
	}

	
	public String getPathwayInteraction(String pathwayname,String organism) {

		List<String> typeIDList = new ArrayList<String>();
		typeIDList.add("");
		typeIDList.add("");
		
		List<String> typeNameList = new ArrayList<String>();
		typeNameList.add("gene_A");
		typeNameList.add("gene_B");
		
		String sqlQuery = "select gene_A,gene_B from gene_pairs where gene_pairs.pathway_name=\""
			+ pathwayname + "\"";
		return getResult(organism, typeIDList, typeNameList, sqlQuery);
	}
	
	public String getPathwayDiff(String pathwayname_A,String pathwayname_B,String organism) {
		
		String gene_A = this.getPathwayGene(pathwayname_A, organism);
		
		String gene_B = this.getPathwayGene(pathwayname_B, organism);
		
		String result = findDiff(gene_A.split("\n"), gene_B.split("\n"));
		
		result += findDiff(gene_B.split("\n"), gene_A.split("\n"));
		
		return result;
	}
	
	public static String findDiff(String[] aArray, String[] bArray){
		
		String diff = "";
		
		for(String a:aArray){
			boolean store = true;
			for(String b:bArray){
				if(a.equalsIgnoreCase(b)){
					store = false;
					break;
				}
			}
			if(store){
				diff += a + "\n";
			}
		}
		return diff;
	}

	public String getIntPathGenes (String organism){
		String typeName = "gene_names";
		String sqlQuery = "select DISTINCT gene_names from pathway_genes";
		return getResult(organism, typeName, sqlQuery);		
	}	
	
	public String getIntPathGenePairs (String organism){
		
		List<String> typeIDList = new ArrayList<String>();
		typeIDList.add("");
		typeIDList.add("");		
		
		List<String> typeNameList = new ArrayList<String>();
		typeNameList.add("gene_A");
		typeNameList.add("gene_B");
		String sqlQuery = "select DISTINCT gene_A,gene_B from gene_pairs";
		
		return getResult(organism, typeIDList, typeNameList, sqlQuery);		
	}
	
	public String getIntPathPathways (String organism){
		String typeName = "pathway_name";
		String sqlQuery = "select DISTINCT pathway_name from pathway_genes";
		return getResult(organism, typeName, sqlQuery);		
	}	
	
	public String getResult(String organism, String typeName, String sqlQuery){
		
		List<String> typeIDList = new ArrayList<String>();
		typeIDList.add("");
		
		List<String> typeNameList = new ArrayList<String>();
		typeNameList.add(typeName);
		
		return getResult(organism, typeIDList, typeNameList, sqlQuery);
	}
	


	public String getResult(String organism, List<String> typeIDList, List<String> typeNameList, String sqlQuery){
		
		StringBuilder result = new StringBuilder("");
		Connection con = null;
		Statement smtm = null;
		ResultSet rs = null;
	
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
			
			rs = smtm.executeQuery(sqlQuery);
			
			while (rs != null && rs.next()) {
				for(int i = 0; i < typeNameList.size(); i++){
					result.append(typeIDList.get(i));
					result.append(rs.getString(typeNameList.get(i)));
					if(i + 1 == typeNameList.size())
						result.append("\n");
					else
						result.append("\t");
				}
			}
			
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
			} catch (Exception ee) {}
		}
	
		return result.toString();
	   }
	}
