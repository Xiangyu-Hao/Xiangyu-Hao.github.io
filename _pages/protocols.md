---
layout: archive
title: "Protocols"
permalink: /protocols/
author_profile: true
---
<script src="https://cdn.jsdelivr.net/npm/js-sha256@0.9.0/src/sha256.min.js"></script>
<script>
    function checkPassword() {
        const obfuscatedHashParts = [
            'e31003f5', '98179cf0', '3921ecde', '0dd03332', 
            '213197ef', '463e77c8', '6e9c9c15', '936c53e8'
        ];
        
        const getCorrectHash = () => {
            return obfuscatedHashParts[7] + obfuscatedHashParts[2] + obfuscatedHashParts[3] + 
                   obfuscatedHashParts[5] + obfuscatedHashParts[0] + obfuscatedHashParts[6] + 
                   obfuscatedHashParts[1] + obfuscatedHashParts[4];
        };
        const correctHash = getCorrectHash();
        
        const userPassword = prompt("请输入访问密码：");
        if (userPassword === null || userPassword === "") {
            return;
        }
        const userHash = sha256(userPassword);
        
        if (userHash === correctHash) {
            window.location.href = "../protocols/Phylogenomics.html";
        } else {
            alert("密码错误，请重试！");
        }
    }
</script>
# <a onclick="checkPassword()">Phylogenomics</a>
# <a href="../protocols/transcriptome.html" target="_blank">Transcriptome</a>
# <a href="../protocols/Molecular Convergence (conv_cal).html" target="_blank">Molecular Convergence (conv_cal)</a>
# <a href="../protocols/Relaxed Selection Test.html" target="_blank">Relaxed Selection Test</a>
