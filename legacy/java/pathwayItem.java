/* Title: IntPath Utilities for Hypergeometric analysis
 * Function: Program for IntPath Hypergeometric utilities
 * 
 * Author: Hufeng Zhou
 * Time: August 9th 2021
 * Version: v2
 */

package utils;

import  java.util.*;
import  java.io.*;

public class pathwayItem {
	
	String name;
	double score;
	int n=0;
	int k=0;
	int m=0;
	int N=0;
	
	public pathwayItem(String name, double score, int n, int k, int m,int N) {
		//super();
		this.name = name;
		this.score = score;
		this.n = n;
		this.k = k;
		this.m = m;
		this.N = N;
	}
	public pathwayItem(String name, double score) {
		super();
		this.name = name;
		this.score = score;
	}
	public pathwayItem() {}
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public double getScore() {
		return score;
	}
	public void setScore(double score) {
		this.score = score;
	}
	public int getn() {
		return n;
	}
	public void setn(int n) {
		this.n = n;
	}
	public int getK() {
		return k;
	}
	public void setK(int k) {
		this.k = k;
	}
	public int getM() {
		return m;
	}
	public void setM(int m) {
		this.m = m;
	}
	public int getN() {
		return N;
	}
	public void setN(int N) {
		this.N = N;
	}
	public void set(String name, double score, int n, int k, int m,int N) 
	{
		//super();
		this.name = name;
		this.score = score;
		this.n = n;
		this.k = k;
		this.m = m;
		this.N = N;
	}
	

}
