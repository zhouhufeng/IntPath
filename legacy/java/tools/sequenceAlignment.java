package tools;

/**
 * This program is used to do sequence alignament of the two strings to see
 * The similarity between pathway names among three different pathway databases.
 * The good thing is when I compare the strings there is no complex matrix.
 * If two characters match, score is 5.
 * If two characters mismatched, score is -2.
 * If two characters match to space(insertion or deletion), score is -1.
 * @author willam
 */


public class sequenceAlignment {
    
//    static String a = "cytochrome C oxidasation";
//    static String b = "miR-1 in cardiac development";
//    static String b = "NADH to cytochrome <i>bd</i> oxidase electron transfer";
    //public static void main(String args[])throws Exception{
    public double thres;
    public int alsc;

    public sequenceAlignment(String a, String b){

        alsc = calcAlignmentScore(a,b);
        
        thres = (double)2*alsc/(a.length()+b.length());

        //System.out.println(aligscore+"\t"+thres);
        
    }

    public int calcAlignmentScore(String a, String b){

        int alignScore = 1;
        a =" "+a;
        b =" "+b;

        int[][] score = new int[a.length()+1][b.length()+1];
        int scoreindel   = 0;
        int scorem  = 0;
        //int scoremism = 0;
        //System.out.println(a.length()+"\t"+b.length());

        for(int i = 0; i< a.length();i++){
            for(int j = 0; j<b.length();j++){

                String cha = a.substring(i, i+1);
                String chb = b.substring(j, j+1);

                scoreindel = 0;

                if(cha.equalsIgnoreCase(chb)){
                    scorem = 1;
                }else{
                    scorem = -1;
                }


                if(i==0 && j==0){
                    
                    score[i][j]=0;

                }else if(i==0){

                    score[i][j]=score[0][j-1]+scoreindel;

                }else if(j==0){

                    score[i][j]=score[i-1][0]+scoreindel;

                }else if(i>0&&j>0){

                    int x =score[i][j-1]+scoreindel;
                    int y =score[i-1][j]+scoreindel;
                    int z =score[i-1][j-1]+scorem;

                    int m = -100;
                    if(x>m) m=x;
                    if(y>m) m=y;
                    if(z>m) m=z;
                    score[i][j]=m;
                    //System.out.print(m+"\t");
                }               
            }
            //System.out.print("\n");
        }
        //System.out.println(score[a.length()-1][b.length()-1]);

        alignScore = score[a.length()-1][b.length()-1];
        //System.out.println(alignScore);
        
        return alignScore;

    }
}
