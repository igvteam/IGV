package org.broad.igv.drs;

import org.broad.igv.oauth.OAuthUtils;
import org.broad.igv.util.HttpUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Map;

// https://drshub.dsde-prod.broadinstitute.org/
public class TerraResolver {


    public static Object resolve(String drsURL) {

        try {
            OAuthUtils.getInstance().getGoogleProvider().checkLogin();

            String accessToken = OAuthUtils.getInstance().getGoogleProvider().getAccessToken();

            System.out.println(accessToken);

            Object response = doPost("https://drshub.dsde-prod.broadinstitute.org/", accessToken);

            return response;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static Object doPost(String drsURL, String accessToken) throws IOException {
        URL url = new URL(drsURL);
        HttpURLConnection con = (HttpURLConnection) url.openConnection();
        con.setRequestMethod("POST");
        con.setRequestProperty("Accept", "*/*");
        con.setRequestProperty("Content-Type", "application/json");
        con.setRequestProperty("Authorization", "Bearer " + accessToken);
        con.setDoOutput(true);

        String body = """
                 {
                 "url": """ +   drsURL + """ 
                "fields": ["bucket", "contentType", "fileName", "gsUri", "hashes", "localizationPath", "name", "size", "timeCreated", "timeUpdated"]
                }       
                """;

        try(OutputStream os = con.getOutputStream()) {
            byte[] input = body.getBytes("utf-8");
            os.write(input, 0, input.length);
        }

        try(BufferedReader br = new BufferedReader(
                new InputStreamReader(con.getInputStream(), "utf-8"))) {
            StringBuilder response = new StringBuilder();
            String responseLine = null;
            while ((responseLine = br.readLine()) != null) {
                response.append(responseLine.trim());
            }
            System.out.println(response.toString());
        }

        return null;
    }

}


//curl -X 'POST' \
//        'https://drshub.dsde-prod.broadinstitute.org/api/v4/drs/resolve' \
//        -H 'accept: */*' \
//        -H 'Authorization: Bearer <token for Terra user>' \
//        -H 'Content-Type: application/json' \
//        -d '{
//        "url": "drs://dg.4503:dg.4503/c64fded3-459b-4fca-82bc-35d154e9aa91",
//        "fields": ["bucket", "contentType", "fileName", "gsUri", "hashes", "localizationPath", "name", "size", "timeCreated", "timeUpdated"]
//        }'


// Compact IDS - not terra
//   drs://dg.4503:dg.4503/c64fded3-459b-4fca-82bc-35d154e9aa91