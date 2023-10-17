#include "wbspeck.h"

void Key_expand32(u16 key[4], u16 roundkey[22])
{
    u16 i, b = key[3];
    u16 a[SPECK_KEY_LEN - 1];

    for(i = 0; i < (SPECK_KEY_LEN - 1); i++)
    {
        a[i] = key[i + 1];
    }
    roundkey[0] = b;  
    for(i = 0; i < SPECK_ROUNDS - 1; i++)
    {
        R32(a[i % (SPECK_KEY_LEN - 1)], b, i);
        roundkey[i + 1] = b;
    }
}
u16 ROR32(u16 x, int r)
{
    return (x >> r) | (x << ((sizeof(u16) * 8) - r));
}
u16 ROL32(u16 x, int r)
{
    return (x << r) | (x >> ((sizeof(u16) * 8) - r));
}

void CP_DCA32()
{
    u16 x0, x1, y0, y1, z0, z1;
    int k, x, j, s, ts, l, b;
    u32 state;
    u16 key[4] = {0x1918, 0x1110, 0x0908, 0x0100};
    u16 roundkey[SPECK_ROUNDS];
    Key_expand32(key, roundkey);
    FILE *fp = NULL;
    ////
    int f = 0, g = 0;
    static const int ktb16 = 0x0000;
    static const int  kte16 = 0xffff;
    int key_count = 0;

    Aff32 Aff, Aff_inv;
    genaffinepairM32(&Aff, &Aff_inv);
    u32 map[pt]; 
    
    u16 key_guess[10000] = {0};
    int k_count = 0;
    double k_max = 0.0;
    int knum = 0;
    for(k = ktb16; k <= kte16; k++) // key
    {
        for(x = 0; x < pt; x++) // adaptive inputs
        {
            x0 = x;
            x1 = 0x00ff;
            
            y1 = ROR32(ROL32(x0 - x1, 7) ^ x1, 2);
            y0 = ROL32((ROL32(x0 - x1, 7) ^ k) - y1, 7);
            
            z1 = ROL32(y1, 2) ^ (ROR32(y0, 7) + y1) ^ roundkey[0];
            z0 = ROR32((ROR32(y0, 7) + y1) ^ roundkey[0], 7) + z1;

            state = (z0 << 16) | z1;
            map[x] = state;
            map[x] = affineU32(Aff, state); // affine encoding
        }
        double k_score = 0.0;
        double L_score[pt] = {0.0};
        int L_count[pt] = {0};
        int l_count = 0;
        double l_max = 0.0;
        for(l = 1; l < pt; l++)
        {  
            double score = 0.0;
            double j_max = 0.0;
            int j_count = 0;
            for(j = 0; j < tb32; j++) // samples of traces
            {
                int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                for(x = 0; x < pt; x++)
                { 
                    if(map[x] & idM32[j]) ts = 1; // traces
                    else ts = 0;
                    ////// correlation
                    f = ts;
                    if(f) Nf1++;
                    else Nf0++;
                    
                    g = xor[l & x];
                    if(g) Ng1++;
                    else Ng0++;
                    
                    if(f == 1 && g == 1) N11++;
                    else if(f == 0 & g == 0) N00++;
                    else if(f == 1 & g == 0) N10++;
                    else N01++;
                }
                if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                else score = 0.0;
                if(score > j_max)
                {
                    j_count = 0;
                    j_max = score;
                    j_count++;
                }
                else if(fabs(score - j_max) < EPS)
                {
                    j_count++;
                }
            }
            L_score[l] = j_max;
            L_count[l] = j_count;
            if(L_score[l] > l_max) l_max = L_score[l];
        }
        k_score = l_max;
        if(fabs(k_score - k_max) < EPS)
        {
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            if(l_count > k_count)
            {
                k_count = l_count;
                knum = 0;
                key_guess[knum] = k;
                knum++;
                printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                fp = fopen("Result_SPECK32.txt", "a");
                fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                fclose(fp);
            }
            else if(l_count == k_count)
            {
                key_guess[knum] = k;
                knum++;
                printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                fp = fopen("Result_SPECK32.txt", "a");
                fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                fclose(fp);
            }
        }
        else if(k_score > k_max) 
        {
            k_max = k_score; 
            knum = 0;
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            k_count = l_count;
            key_guess[knum] = k;
            knum++;
            printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            fp = fopen("Result_SPECK32.txt", "a");
            fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            fclose(fp);
        }
    }
    printf("------\n"); 
    fp = fopen("Result_SPECK32.txt", "a");
    fprintf(fp, "------\n"); 
    fclose(fp); 
    for(k = 0; k < knum; k++)
    {
        printf("simulation WB-SPECK32 recoverd key: %.2x, correlation: %f, encoding count: %d\n", key_guess[k], k_max, k_count);
        fp = fopen("Result_SPECK32.txt", "a");
        fprintf(fp, "simulation WB-SPECK32 recoverd key: %.2x, correlation: %f, encoding count: %d\n", key_guess[k], k_max, k_count);
        fclose(fp);
        key_count++;
    }
    printf("------\n"); 
    fp = fopen("Result_SPECK32.txt", "a");
    fprintf(fp, "------\n"); 
    fclose(fp); 

    printf("simulation WB-SPECK32 recoverd key count: %d\n", key_count);
    fp = fopen("Result_SPECK32.txt", "a");
    fprintf(fp, "simulation WB-SPECK32 recoverd key count: %d\n", key_count);
    fclose(fp); 
}
