package converter;

import java.util.Arrays;

import static java.lang.Math.*;

public class CoordinateConverter {

    private static final double LAMBDA0 = 0.33246029531;
    private static final double E = 0.0818205679407;
    private static final double K = 1.003110007693;
    private static final double N = 1.000719704936;
    private static final double FI0 = 0.822050077;
    private static final double R = 6379743.001;
    private static final double M0 = 0.99993;

    private double[] toGaussSpherical(double x, double y) {
        double xInRadian = toRadians(x);
        double yInRadian = toRadians(y);
        double u = ((PI / 4) + (xInRadian / 2));
        double p =  pow(((1 - (E * sin(xInRadian))) / (1 + (E * sin(xInRadian)))), (N * E) / 2);
        double k = atan(K * pow(tan(u), N) * p);
        double fi = 2 * (k  - (PI / 4));
        double lambda = N * (yInRadian - LAMBDA0);
        return new double[]{fi, lambda};
    }

    private double[] toAuxSpherical(double x, double y) {
        double fi = toGaussSpherical(x, y)[0];
        double lambda = toGaussSpherical(x, y)[1];
        double fiComma = asin((cos(FI0) * sin(fi)) - (sin(FI0) * cos(fi) * cos(lambda)));
        double lambdaComma = asin((cos(fi) * sin(fi)) / cos(fiComma));
        return new double[]{fiComma, lambdaComma};
    }

    private double[] toAuxFlat(double x, double y) {
        double fiComma = toAuxSpherical(x, y)[0];
        double lambdaComma = toAuxSpherical(x, y)[1];
        double xT = R * M0 * log(tan((PI / 4) + (fiComma / 2)));
        double yT = R * M0 * lambdaComma;
        return new double[]{xT, yT};
    }

    public double[] wgs84ToEov(double x, double y){
       double[] auxFlatCoordinates = toAuxFlat(x, y);
       return new double[]{
               auxFlatCoordinates[0] + 200000,
               auxFlatCoordinates[1] + 650000
       };
    }

    public static void main(String[] args) {
        CoordinateConverter con = new CoordinateConverter();

        System.out.println(Arrays.toString(con.wgs84ToEov(46.264748, 20.157790)));
    }
}
