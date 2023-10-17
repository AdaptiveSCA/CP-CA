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

void CP_ADCA32()
{
    u16 x0, x1, y0, y1, z0, z1;
    int k, x, j, s, ts, l, b, bit;
    u32 state;
    u16 key[4] = {0x1918, 0x1110, 0x0908, 0x0100};
    u16 roundkey[SPECK_ROUNDS];
    Key_expand32(key, roundkey);
    ////
    int f = 0, g = 0;
    static const int ktb16 = 0x0000;
    static const int  kte16 = 0xffff;
    int key_count = 0;

    Aff32 Aff, Aff_inv;
    genaffinepairM32(&Aff, &Aff_inv);
    u32 map[pt]; 
    FILE *fp = NULL;

    uint8_t Htrace[pt][9];
    uint8_t temp;
    uint8_t trail[1280][3];// Gaussian trail
    for(x = 0; x < pt; x++)
    {
        Htrace[x][0] = 1;
        for(bit = 1; bit < 9; bit++)
        {
            if(x & idM8[bit - 1]) Htrace[x][bit] = 1;
            else Htrace[x][bit] = 0;
        }
    }
    int Gauss_time = 0; // 1277
    //  Gaussian Elimination
    for(int i = 0; i < 9; i++)
    {
        if(Htrace[i][i])
        {
            for(int j = i + 1; j < pt; j++)
            {
                if(Htrace[j][i])
                { 
                    for(int r = 0; r < 9; r++)
                    {
                        Htrace[j][r] ^= Htrace[i][r];
                    }
                    trail[Gauss_time][0] = 1; //addition
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;
                }
            }
        }
        else
        {
            int flag = 0;
            for(int j = i + 1; j < pt; j++)
            {
                if(Htrace[j][i])
                {
                    for(int r = 0; r < 9; r++)
                    {
                        temp = Htrace[i][r];
                        Htrace[i][r] = Htrace[j][r];
                        Htrace[j][r] = temp;
                    }
                    trail[Gauss_time][0] = 0; //swap
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;

                    flag = 1;
                    break;
                }
            }
            if(flag)
            {
                for(int j = i + 1; j < pt; j++)
                {
                    if(Htrace[j][i])
                    { 
                        for(int r = 0; r < 9; r++)
                        {
                            Htrace[j][r] ^= Htrace[i][r];
                        }
                        trail[Gauss_time][0] = 1; //addition
                        trail[Gauss_time][1] = j;
                        trail[Gauss_time][2] = i;
                        Gauss_time++;
                    }
                }
            }
        }
    }
    // printf("Gauss_time: %d\n", Gauss_time);
    // return ;

    // the rank of A
    int rA = pt;
    for(int i = pt - 1; i >= 0; i--)
    {
        int allzero = 1;
        for(int j = 0; j < 9; j++)
        {
            if(Htrace[i][j]) allzero = 0;
        }
        if(allzero) rA--;
        else break;
    }
    // printf("rA: %d\n", rA);
    // return ;

    u16 key_guess[10000] = {0};
    int k_max = 0;
    int knum = 0;
    uint8_t vector[pt];
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
        int k_count = 0;
        for(j = 0; j < tb32; j++) // samples of traces
        {
            for(x = 0; x < pt; x++)
            { 
                if(map[x] & idM32[j]) vector[x] = 1; // traces
                else vector[x] = 0;
            }
            //  Gaussian Elimination
            for(int i = 0; i < Gauss_time; i++)
            {
                if(trail[i][0]) // addition
                {
                    vector[trail[i][1]] ^= vector[trail[i][2]];
                }
                else // swap
                {
                    temp = vector[trail[i][2]];
                    vector[trail[i][2]] = vector[trail[i][1]];
                    vector[trail[i][1]] = temp;
                }
            }
            
            // Gauss Over
            int rAb = pt;
            for(int i = pt - 1; i >= 0; i--)
            {
                int allzero = 1;
                if(vector[i]) allzero = 0;
                if(allzero) rAb--;
                else break;
            }
            if(rA >= rAb) // has a solusion
            {
                k_count++;
            }
        }
        if(k_count && (k_count == k_max))
        {
            key_guess[knum] = k;
            knum++;
            fp = fopen("Result_SPECK32.txt", "a");
            fprintf(fp, "key guess: %.2x, encoding count: %d\n", k, k_count);
            printf("key guess: %.2x, encoding count: %d\n", k, k_count);
            fclose(fp);
        }
        else if(k_count && (k_count > k_max))
        {
            k_max = k_count; 
            knum = 0;
            key_guess[knum] = k;
            knum++;
            fp = fopen("Result_SPECK32.txt", "a");
            fprintf(fp, "key guess: %.2x, encoding count: %d\n", k, k_count);
            printf("key guess: %.2x, encoding count: %d\n", k, k_count);
            fclose(fp);
        }
    }
    fp = fopen("Result_SPECK32.txt", "a");
    fprintf(fp, "------\n");
    printf("------\n");
    fclose(fp);  
    for(k = 0; k < knum; k++)
    {
        fp = fopen("Result_SPECK32.txt", "a");
        fprintf(fp, "simulation WB-SPECK32 recoverd key: %.2x, encoding count: %d\n", key_guess[k], k_max);
        printf("simulation WB-SPECK32 recoverd key: %.2x, encoding count: %d\n", key_guess[k], k_max);
        fclose(fp); 
        key_count++;
    }
    fp = fopen("Result_SPECK32.txt", "a");
    fprintf(fp, "------\n");
    printf("------\n");
    fclose(fp);  

    fp = fopen("Result_SPECK32.txt", "a");
    fprintf(fp, "simulation WB-SPECK32 recoverd key count: %d\n", key_count);
    printf("simulation WB-SPECK32 recoverd key count: %d\n", key_count);
    fclose(fp); 
}
