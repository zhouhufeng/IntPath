/* Title: IntPath Pathway Analysis
 * Function: Program for IntPath Pathway Enrichment Analysis 
 * 
 * Author: Hufeng Zhou
 * Time: August 9th 2021
 * Version: v2
 */

package utils;
import java.io.*;
import java.math.*;
import java.util.*;
import java.sql.*;

public class GeneListsAnalyse {

	public String organism; 
	
	public int[][] dist;		

	public GeneListsAnalyse() {

	}	
	
	public void setOrganism(String organism) {
		this.organism = organism;
	}

	public Vector<String> setGeneList(String text, Vector<String> geneList) {
		
		String[] line_temp = text.split("\r\n");			

		for (int i = 0; i < line_temp.length; i++)
			geneList.add(line_temp[i]);
		
		return geneList;
	}

	public void setGeneList_file(String path, Vector<String> geneList) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(path));

		String readline;

		while (br.ready()) {
			readline = br.readLine();
			String[] line_temp = readline.split("\t");
			geneList.add(line_temp[0]);
		}
		
	}

	public int getsize(Vector<String> geneList ) {
		// String[] line = this.text.split("\n");
		return geneList.size();
	}

	public static BigInteger choose(int n, int m) {

		if (m > n / 2)
			return choose(n, n - m);

		BigInteger r = new BigInteger("1");

		for (int i = 0; i < m; i++)
			r = r.multiply(new BigInteger("" + (n - i))).divide(new BigInteger("" + (i + 1)));
		return r;
	}

	public static double hyperGeometric(int N, int n, int m, int k) {
		return 1.0 / (choose(N, n).divide(choose(m, k)).divide(choose(N - m, n- k))).doubleValue();
	}

	public static void mergesortpathwayItem_x(Vector<pathwayItem> sequence,	int m, int n) {
		
		if (m == n)
			return;
		int k = (m + n) / 2;

		if (sequence.size() == 0) {
			System.out.println("Empty");
		} else {
			System.out.println("size: ");
			System.out.println(sequence.size());
			mergesortpathwayItem_x(sequence, m, k);
			mergesortpathwayItem_x(sequence, k + 1, n);
			mergepathwayItem_x(sequence, m, k, n);
		}
	}

	public static void mergepathwayItem_x(Vector<pathwayItem> sequence, int m, int k, int n) {

		Vector<pathwayItem> temp = new Vector<pathwayItem>();

		for (int i = 0; i < n - m + 1; i++) {
			pathwayItem temp1 = new pathwayItem();
			temp.add(temp1);
		}

		int p = m;
		int q = k + 1;
		int r = 0;

		while (p <= k && q <= n) {

			if (sequence.get(p).getScore() <= sequence.get(q).getScore()) {
				temp.set(r, sequence.get(p));
				p++;
				r++;
			} else {
				temp.set(r, sequence.get(q));
				q++;
				r++;
			}
		}

		while (p <= k) {
			temp.set(r, sequence.get(p));
			p++;
			r++;

		}

		while (q <= n) {
			temp.set(r, sequence.get(q));
			q++;
			r++;
		}

		for (int i = m; i <= n; i++) {
			sequence.set(i, temp.get(i - m));
		}

	}

	public void sort(Vector<pathwayItem> pathways) {
		mergesortpathwayItem_x(pathways, 0, pathways.size() - 1);
	}
	
	public Vector<pathwayItem> findPaths (double threshold,Vector<String> geneList, Vector<pathwayItem> pathways) throws IOException {

		Vector<pathwayItem> pathway_temp = new Vector<pathwayItem>();

		Vector<String> gene_list = new Vector<String>(); 

		for (int i = 0; i < geneList.size(); i++) {

			if (!gene_list.contains(geneList.get(i))) {

				gene_list.add(geneList.get(i));

			}
		}

		int N = SpeciesInformation.getGenesCount(organism);

		int n = gene_list.size();

		Vector<String> pathway_trace = new Vector<String>();

		Connection con = null;

		Statement smtm = null;

		ResultSet rs = null;

		try {

			String databaseURI = SpeciesInformation.getDatabaseURI(organism);

			Class.forName("com.mysql.jdbc.Driver");

			con = DriverManager.getConnection(databaseURI);

			smtm = con.createStatement();

			for (int i = 0; i < gene_list.size(); i++) {

				String query = "select pathway_name from pathway_genes where  pathway_genes.gene_names =\""	+ gene_list.get(i) + "\"";

				rs = smtm.executeQuery(query);

				while (rs != null && rs.next()) {

					int m = 0, k = 0;

					String pathway_names = rs.getString("pathway_name");
					// System.out.println(pathway_id);

					int count = 0;

					int signal = 1;

					while (count < pathway_trace.size()) {

						if (pathway_names.equals(pathway_trace.get(count))) {
							signal = 0;
							break;
						} else
							count++;
					}
					
					if (signal == 1) {
						
						pathway_trace.add(pathway_names);						

						Statement smtm2 = null;

						ResultSet rs2 = null;

						String query2 = "select gene_names from pathway_genes where pathway_genes.pathway_name=\""+ pathway_names + "\"";

						smtm2 = con.createStatement();

						rs2 = smtm2.executeQuery(query2);

						while (rs2 != null && rs2.next()) {

							m++;

							String gene_name = rs2.getString("gene_names");

							if (gene_list.contains(gene_name)) {
								k++;
								//break;//hidden bug! left by jingjing
							}
						}

						double sum = 0;

						for (int j = k; j <= n; j++) {
							double hyper_test = hyperGeometric(N, n, m, j);
							sum = sum + hyper_test;
						}

						pathwayItem temp = new pathwayItem(pathway_names, sum, n, k, m, N);

						if (sum <= threshold)
							pathway_temp.add(temp);

					}
				}
			}
			
		}

		catch (ClassNotFoundException e) {
			e.printStackTrace();
		}

		catch (SQLException ex) {

			System.out.println(ex.getMessage());
			System.out.println(ex.getSQLState());
			System.out.println(ex.getErrorCode());

		}

		finally {

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
		return pathway_temp;
	}
	
	public int[][] getDist() {
		return dist;
	}

	public void setDist(int[][] dist) {
		this.dist = dist;
	}
			
	public void distance(String rootPath, Vector<pathwayItem> pathway1,Vector<pathwayItem> pathway2) throws Exception
	{
		FloydWarshall f = new FloydWarshall();
		System.out.print(rootPath);
		
		String DISTdir = SpeciesInformation.getOrgDIST(organism);

		f.loadModel(rootPath+DISTdir); 

		dist = new int[pathway1.size()][pathway2.size()];
		
		Connection con = null; 
		
        Statement smtm = null;
        
        ResultSet rs   = null;
        
		try
        {
		 String databaseURI = SpeciesInformation.getDatabaseURI(organism);

		 Class.forName("com.mysql.jdbc.Driver");

		 con = DriverManager.getConnection(databaseURI);
			
		 smtm=con.createStatement();
          
          for(int i=0;i<pathway1.size();i++)
          {
        	  String query1= "select gene_names from pathway_genes where pathway_genes.pathway_name=\""+pathway1.get(i).getName()+"\""; 
        	  
        	  rs=smtm.executeQuery(query1);
        	  
        	  Vector<String> h1=new Vector<String>();
        	  
        	  while(rs.next())
        	  {
        		  String gene_name=rs.getString("gene_names");
        		  
        		  h1.add(gene_name);
        		  
        		  //System.out.println("gene name h2 debug:"+gene_name);
        	  }
        	  
        	  Statement smtm1=null;
              ResultSet rs1=null;
              smtm1=con.createStatement();
              
              for(int j=0;j<pathway2.size();j++)
              {
            	  String query2= "select gene_names from pathway_genes where pathway_genes.pathway_name=\""+pathway2.get(j).getName()+"\"";
            	  
            	  rs1=smtm1.executeQuery(query2);
            	  
            	  Vector<String> h2=new Vector<String>();
            	  
            	  while(rs1.next())
            	  {
            		  String gene_name=rs1.getString("gene_names");
            		  
            		  h2.add(gene_name);
            		  
            		  //System.out.println("gene name h2 debug:"+gene_name);
            	  }
            	  
            	  int d=Integer.MAX_VALUE;
            	  
            	  double sum=0;
            	  
            	  int num=0;
            	  
            	  for(int p=0;p<h1.size();p++)
            	  {
            		  for(int q=0;q<h2.size();q++)
            		  {
            			  int temp=f.getDistance(h1.get(p), h2.get(q));
            			  
            			  if(temp<Integer.MAX_VALUE)
            			  {
            			      sum=sum+temp;
            			      num++;
            			  }
            			  if(temp<d)
            				  d=temp;
            		  }
            	  }
            	  if(num>0)
            	       sum=sum/num;
            	  
            	  if(d<Integer.MAX_VALUE)
            	  {
            	       dist[i][j]=d;
            	  }
            	  else
            	  {
            		   dist[i][j]=-1;
            	  }
              }  
              if (rs1!=null)
              {
                 rs1.close();
                 rs1=null;
              }
          
              if(smtm1!=null)
              {
                 smtm1.close();
                 smtm1=null;
           
              } 
          }
           
        }
        
        catch(ClassNotFoundException e) 
        {
         e.printStackTrace();
        }
        
        catch(SQLException ex)
        {
         System.out.println(ex.getMessage());
         System.out.println(ex.getSQLState());
         System.out.println(ex.getErrorCode());        
        }
           try
           {
              if (rs!=null)
             {
                rs.close();
                rs=null;
             }
         
             if(smtm!=null)
             {
                smtm.close();
                smtm=null;
          
             }
         
             if(con!=null)
             {
               con.close();
               con=null;
             }
             
           }
           catch(Exception ee)
           {
         
           }
		
	}

}
